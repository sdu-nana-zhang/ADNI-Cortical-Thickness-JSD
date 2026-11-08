source("src/project.R")
setwd(stage_dir("edges"))
suppressWarnings(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(REPO_ROOT,'src/json_io.R'))
read<-function(f)read.csv(f,fileEncoding='UTF-8-BOM',check.names=FALSE,stringsAsFactors=FALSE)
write<-function(x,f)write.csv(x,f,row.names=FALSE,fileEncoding='UTF-8',na='')
cfg<-fromJSON('protocol/frozen_config.json');G<-cfg$n_patients;N<-cfg$n_MRI;ed<-read('inputs/edge_dictionary.csv');rois<-read('inputs/ROI_dictionary.csv')
files<-sprintf('models/edge_%04d.rds',1:1891);stopifnot(all(file.exists(files)))
recs<-lapply(files,readRDS);stopifnot(all(vapply(seq_along(recs),function(k)recs[[k]]$EDGE_ID==k&&recs[[k]]$protocol_sha256==cfg$protocol_sha256&&recs[[k]]$implementation_id==cfg$implementation_id,logical(1))))
fail<-list();nr<-0
for(k in 1:1891)for(m in c('L0','L1'))if(!isTRUE(recs[[k]]$models[[m]]$valid)){nr<-nr+1;fail[[nr]]<-data.frame(EDGE_ID=k,MODEL=m,ERROR=recs[[k]]$models[[m]]$error)}
if(nr>0){write(do.call(rbind,fail),'audit/incomplete_models.csv');stop(paste(nr,'edge models do not have valid fits; do not silently reduce the inference family'))}
res<-list();coefs<-list();vars<-list();diagnostic<-list();ii<-0
signature<-'latest_ADNI_770_model'
for(k in 1:1891){
 stopifnot(identical(recs[[k]]$input_signature,signature))
 for(m in c('L0','L1')){
  a<-recs[[k]]$models[[m]];ii<-ii+1;res[[ii]]<-a$stats
  coefs[[ii]]<-data.frame(EDGE_ID=k,MODEL=m,TERM=names(a$beta),BETA=as.numeric(a$beta),CR1_SE=sqrt(diag(a$CR1_cov)),MODEL_SE=sqrt(diag(a$model_cov)),row.names=NULL)
  vars[[ii]]<-data.frame(EDGE_ID=k,MODEL=m,SIGMA2=a$sigma2,D00=a$D00,D11=a$D11,RANDOM_TYPE=a$random_type,N_PATIENTS=G,N_MRI=N,P=length(a$beta),DF=G-1,CR1_FACTOR=G/(G-1)*(N-1)/(N-length(a$beta)),row.names=NULL)
  diagnostic[[ii]]<-data.frame(EDGE_ID=k,MODEL=m,RANDOM_TYPE=a$random_type,OPTIMIZER=a$optimizer,ATTEMPTS=length(a$attempts),NEAR_BOUNDARY=a$near_boundary,GLS_BETA_MAX_ERROR=a$GLS_beta_error,GLS_COV_MAX_ERROR=a$GLS_cov_error,LOGLIK=a$logLik,row.names=NULL)
 }
}
res<-do.call(rbind,res);res$Q_BH1891<-NA_real_
for(m in c('L0','L1'))for(cc in c('CN','AD','AD-CN')){
 i<-which(res$MODEL==m&res$CONTRAST==cc);stopifnot(length(i)==1891,all(is.finite(res$P[i])));res$Q_BH1891[i]<-p.adjust(res$P[i],method='BH',n=1891)
}
res$FAMILY<-paste(res$MODEL,res$CONTRAST,sep=':');res$PRIMARY_FAMILY<-res$MODEL=='L0'&res$CONTRAST%in%c('AD-CN');res$SUPPORTED_Q05<-res$Q_BH1891<.05
res$DIRECTION<-ifelse(res$TYPE=='GROUP_SLOPE',ifelse(res$ESTIMATE<0,'趋同','分化'),ifelse(res$ESTIMATE<0,'较低年变化率','较高年变化率'))
key<-paste(res$MODEL,res$EDGE_ID,res$CONTRAST,sep='|');res$GROUP1_RATE<-NA_real_;res$GROUP0_RATE<-NA_real_;res$GROUP1_Q<-NA_real_;res$GROUP0_Q<-NA_real_;res$ABSOLUTE_PATTERN<-NA_character_
for(cc in c('AD-CN')){
 i<-which(res$CONTRAST==cc);pair<-strsplit(cc,'-',fixed=TRUE)[[1]];i1<-match(paste(res$MODEL[i],res$EDGE_ID[i],pair[1],sep='|'),key);i0<-match(paste(res$MODEL[i],res$EDGE_ID[i],pair[2],sep='|'),key)
 res$GROUP1_RATE[i]<-res$ESTIMATE[i1];res$GROUP0_RATE[i]<-res$ESTIMATE[i0];res$GROUP1_Q[i]<-res$Q_BH1891[i1];res$GROUP0_Q[i]<-res$Q_BH1891[i0]
 res$ABSOLUTE_PATTERN[i]<-paste0(pair[1],ifelse(res$ESTIMATE[i1]<0,'趋同','分化'),' / ',pair[2],ifelse(res$ESTIMATE[i0]<0,'趋同','分化'))
}
res<-merge(res,ed,by='EDGE_ID',all.x=TRUE,sort=FALSE);res$ROI1_CN<-rois$ROI_CN[match(res$ROI1,rois$ROI)];res$ROI2_CN<-rois$ROI_CN[match(res$ROI2,rois$ROI)];res<-res[order(res$MODEL,res$CONTRAST,res$EDGE_ID),]
write(res,'results/all_edge_effects.csv');write(res[res$MODEL=='L0',],'results/L0_effects.csv');write(res[res$MODEL=='L1',],'results/L1_effects.csv');write(res[res$SUPPORTED_Q05,],'results/significant_effects.csv')
write(do.call(rbind,coefs),'results/fixed_coefficients.csv');write(do.call(rbind,vars),'audit/variance_parameters.csv');diag<-do.call(rbind,diagnostic);write(diag,'audit/model_diagnostics.csv')
for(m in c('L0','L1')){
 P<-length(recs[[1]]$models[[m]]$beta)
 for(kind in c('CR1_cov','model_cov')){
  arr<-array(unlist(lapply(recs,function(a)a$models[[m]][[kind]]),use.names=FALSE),dim=c(P,P,1891));con<-file(paste0('audit/',m,'_',kind,'_float64.bin'),'wb');writeBin(as.double(arr),con,size=8,endian='little');close(con)
 }
}
counts<-list();ci<-0
for(m in c('L0','L1'))for(cc in c('CN','AD','AD-CN'))for(sign in c(-1,1)){
 z<-res[res$MODEL==m&res$CONTRAST==cc,];ci<-ci+1;counts[[ci]]<-data.frame(MODEL=m,CONTRAST=cc,TYPE=z$TYPE[1],DIRECTION=if(z$TYPE[1]=='GROUP_SLOPE')if(sign<0)'趋同' else '分化' else if(sign<0)'较低年变化率' else '较高年变化率',N_ESTIMATED_DIRECTION=sum(sign*z$ESTIMATE>0),N_Q05=sum(z$Q_BH1891<.05&sign*z$ESTIMATE>0),FAMILY_N=1891)
}
counts<-do.call(rbind,counts);write(counts,'results/direction_counts.csv')
for(cc in c('CN','AD','AD-CN')){
 a<-res[res$MODEL=='L0'&res$CONTRAST==cc,];b<-res[res$MODEL=='L1'&res$CONTRAST==cc,];j<-match(a$EDGE_ID,b$EDGE_ID)
 a$L1_ESTIMATE<-b$ESTIMATE[j];a$L1_Q<-b$Q_BH1891[j];a$BOTH_Q05_SAME_DIRECTION<-a$Q_BH1891<.05&a$L1_Q<.05&sign(a$ESTIMATE)==sign(a$L1_ESTIMATE)
 write(a,paste0('results/L0_L1_comparison_',gsub('-','_',cc),'.csv'))
}
write_json(list(n_edges=1891,valid_edge_models=3782,all_families_n1891=TRUE,no_omnibus_gate=TRUE,RI_fallback_models=sum(diag$RANDOM_TYPE=='INTERCEPT_ONLY'),optim_models=sum(diag$OPTIMIZER=='optim'),near_boundary_models=sum(diag$NEAR_BOUNDARY),max_internal_GLS_beta_error=max(diag$GLS_BETA_MAX_ERROR),max_internal_GLS_cov_error=max(diag$GLS_COV_MAX_ERROR)), 'audit/fit_summary.json',auto_unbox=TRUE,pretty=TRUE,digits=17)
sink('audit/R_sessionInfo.txt');sessionInfo();sink();cat('Collected all 3782 valid edge models.\n');print(counts[counts$MODEL=='L0',],row.names=FALSE)
