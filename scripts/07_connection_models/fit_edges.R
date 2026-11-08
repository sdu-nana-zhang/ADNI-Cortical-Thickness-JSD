source("src/project.R")
setwd(stage_dir("edges"))
suppressWarnings(Sys.setlocale('LC_CTYPE','.UTF-8'))
suppressPackageStartupMessages(library(nlme));source(file.path(REPO_ROOT,'src/json_io.R'))
root<-getwd();args<-commandArgs(trailingOnly=TRUE);worker<-if(length(args)>=1)as.integer(args[1]) else 1L;workers<-if(length(args)>=2)as.integer(args[2]) else 1L
cfg<-fromJSON('protocol/frozen_config.json');implementation<-cfg$implementation_id
input_signature<-'latest_ADNI_770_model'
read<-function(f)read.csv(f,fileEncoding='UTF-8-BOM',check.names=FALSE,stringsAsFactors=FALSE)
d<-read('inputs/covariates.csv');ed<-read('inputs/edge_dictionary.csv');N<-nrow(d);E<-nrow(ed)
stopifnot(file.info('inputs/JSD_columnmajor_float64.bin')$size==8*N*E)
con<-file('inputs/JSD_columnmajor_float64.bin','rb');Y<-matrix(readBin(con,'double',n=N*E,size=8,endian='little'),N,E);close(con)
stopifnot(nrow(d)==N,ncol(Y)==E,all(is.finite(Y)),nrow(ed)==E,all(ed$EDGE_ID==seq_len(E)))
d$BASELINE_DX<-factor(d$BASELINE_DX,levels=c('CN','AD'));d$PTID<-factor(d$PTID);d$Y<-0
stopifnot(nlevels(d$PTID)==cfg$n_patients,all(d$BASELINE_DX%in%c('CN','AD')))
id<-d$PTID;n<-as.numeric(rowsum(rep(1,nrow(d)),id));st<-as.numeric(rowsum(d$TIME,id));st2<-as.numeric(rowsum(d$TIME^2,id));G<-length(n);N<-nrow(d);df<-G-1;crit<-qt(.975,df)
forms<-list(L0=Y~TIME*BASELINE_DX+(AGE_C+MALE+EDUC_C)*TIME,L1=Y~TIME*BASELINE_DX+(AGE_C+MALE+EDUC_C)*TIME+CORTEX_BETWEEN+CORTEX_WITHIN+THICKNESS_BETWEEN+THICKNESS_WITHIN)
fixed<-list();pooled<-d[!duplicated(d$PTID),]
for(model in names(forms)){
 X<-model.matrix(delete.response(terms(forms[[model]])),d);P<-ncol(X);S0<-rowsum(X,id);S1<-rowsum(X*d$TIME,id)
 ref<-d[1,,drop=FALSE];ref$AGE_C<-mean(pooled$AGE_C);ref$MALE<-mean(pooled$MALE);ref$EDUC_C<-mean(pooled$EDUC_C);ref$CORTEX_BETWEEN<-mean(pooled$CORTEX_BETWEEN);ref$THICKNESS_BETWEEN<-mean(pooled$THICKNESS_BETWEEN);ref$CORTEX_WITHIN<-0;ref$THICKNESS_WITHIN<-0
 ll<-list()
 for(g in c('CN','AD')){
  r<-ref[rep(1,2),];r$TIME<-c(0,1);r$BASELINE_DX<-factor(rep(g,2),levels=c('CN','AD'));mx<-model.matrix(delete.response(terms(forms[[model]])),r);ll[[g]]<-mx[2,]-mx[1,]
 }
 L<-rbind(CN=ll$CN,AD=ll$AD,'AD-CN'=ll$AD-ll$CN)
 fixed[[model]]<-list(X=X,S0=S0,S1=S1,XTX=crossprod(X),P=P,L=L,CR1=G/(G-1)*(N-1)/(N-P))
 if(worker==1){write.csv(data.frame(PTID=as.character(d$PTID),X,check.names=FALSE),paste0('inputs/design_',model,'.csv'),row.names=FALSE,fileEncoding='UTF-8');write.csv(data.frame(CONTRAST=rownames(L),L,check.names=FALSE),paste0('inputs/contrasts_',model,'.csv'),row.names=FALSE,fileEncoding='UTF-8')}
}
fit_one<-function(k,model){
 dd<-d;dd$Y<-Y[,k];attempts<-list();fit<-NULL;kindchosen<-NULL;optchosen<-NULL
 for(kind in c('DIAGONAL_INTERCEPT_SLOPE','INTERCEPT_ONLY')){
  for(opt in c('nlminb','optim')){
   err<-NULL;ww<-character()
   ans<-tryCatch(withCallingHandlers(lme(fixed=forms[[model]],data=dd,random=if(kind=='DIAGONAL_INTERCEPT_SLOPE')list(PTID=pdDiag(~TIME)) else ~1|PTID,method='REML',na.action=na.fail,control=lmeControl(opt=opt,maxIter=200,msMaxIter=200,niterEM=30,tolerance=1e-6,returnObject=FALSE)),warning=function(w){ww<<-c(ww,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){err<<-conditionMessage(e);NULL})
   attempts[[length(attempts)+1]]<-list(random=kind,optimizer=opt,error=if(is.null(err))'' else err,warnings=ww)
   if(!is.null(ans)&&all(is.finite(fixef(ans)))&&is.finite(ans$sigma)&&!any(grepl('converg|iteration limit',ww,ignore.case=TRUE))){fit<-ans;kindchosen<-kind;optchosen<-opt;break}
  }
  if(!is.null(fit))break
 }
 if(is.null(fit))return(list(valid=FALSE,attempts=attempts,error='No valid converged model'))
 b<-fixef(fit);ob<-fixed[[model]];X<-ob$X;stopifnot(identical(names(b),colnames(X)))
 D<-as.matrix(getVarCov(fit,type='random.effects'));d0<-D[1,1];d1<-if(ncol(D)>1)D[2,2] else 0;s<-fit$sigma^2
 det<-1+(d0*n+d1*st2)/s+d0*d1*(n*st2-st^2)/s^2
 k00<-d0*(1+d1*st2/s)/det;k01<- -d0*d1*st/s/det;k11<-d1*(1+d0*n/s)/det
 A<-ob$XTX/s-(crossprod(ob$S0,ob$S0*k00+ob$S1*k01)+crossprod(ob$S1,ob$S0*k01+ob$S1*k11))/s^2
 bread<-solve(A);dimnames(bread)<-list(names(b),names(b));r<-dd$Y-as.numeric(X%*%b)
 sr<-as.numeric(rowsum(r,id));str<-as.numeric(rowsum(r*dd$TIME,id));U<-rowsum(X*r,id)/s-(ob$S0*(k00*sr+k01*str)+ob$S1*(k01*sr+k11*str))/s^2
 V<-ob$CR1*bread%*%crossprod(U)%*%bread;V<-(V+t(V))/2;dimnames(V)<-list(names(b),names(b))
 sy<-as.numeric(rowsum(dd$Y,id));sty<-as.numeric(rowsum(dd$Y*dd$TIME,id));rhs<-crossprod(X,dd$Y)/s-(crossprod(ob$S0,k00*sy+k01*sty)+crossprod(ob$S1,k01*sy+k11*sty))/s^2
 betaerr<-max(abs(as.numeric(bread%*%rhs)-b));coverr<-max(abs(bread-fit$varFix))
 stopifnot(all(is.finite(V)),all(diag(V)>0),betaerr<1e-7,coverr<1e-7)
 L<-ob$L;est<-as.numeric(L%*%b);se<-sqrt(rowSums((L%*%V)*L));stopifnot(all(is.finite(se)),all(se>0))
 stats<-data.frame(EDGE_ID=k,MODEL=model,CONTRAST=rownames(L),TYPE=c(rep('GROUP_SLOPE',2),'SLOPE_DIFFERENCE'),ESTIMATE=est,CR1_SE=se,CI_LOW=est-crit*se,CI_HIGH=est+crit*se,P=2*pt(-abs(est/se),df),DF=df,N_PATIENTS=G,N_MRI=N,row.names=NULL)
 list(valid=TRUE,beta=b,CR1_cov=V,model_cov=fit$varFix,stats=stats,sigma2=s,D00=d0,D11=d1,random_type=kindchosen,optimizer=optchosen,attempts=attempts,logLik=as.numeric(logLik(fit)),near_boundary=d1/max(s,1e-300)<1e-8,GLS_beta_error=betaerr,GLS_cov_error=coverr)
}
indices<-seq(worker,E,by=workers);done<-0;errors<-0
progress<-function(status,last){write_json(list(worker=worker,workers=workers,status=status,total_edges=length(indices),done_edges=done,error_edges=errors,last_edge=last,utc=format(Sys.time(),tz='UTC',usetz=TRUE)),sprintf('work/progress_%02d.json',worker),auto_unbox=TRUE,pretty=TRUE)}
progress('RUNNING',0)
for(k in indices){
 path<-sprintf('models/edge_%04d.rds',k);old<-if(file.exists(path))readRDS(path) else NULL
 if(!is.null(old)&&identical(old$protocol_sha256,cfg$protocol_sha256)&&identical(old$implementation_id,implementation)&&identical(old$input_signature,input_signature)&&all(vapply(old$models,function(v)isTRUE(v$valid),logical(1)))){done<-done+1;next}
 rec<-list(EDGE_ID=k,FEATURE=ed$FEATURE[k],protocol_sha256=cfg$protocol_sha256,implementation_id=implementation,input_signature=input_signature,models=list())
 for(model in names(forms))rec$models[[model]]<-tryCatch(fit_one(k,model),error=function(e)list(valid=FALSE,error=conditionMessage(e)))
 if(any(!vapply(rec$models,function(v)isTRUE(v$valid),logical(1))))errors<-errors+1
 saveRDS(rec,path,compress='gzip');done<-done+1
 if(done%%10==0||done==length(indices)){progress('RUNNING',k);cat('worker',worker,':',done,'/',length(indices),'edges; errors',errors,'\n');flush.console()}
}
progress(if(errors==0)'COMPLETE' else 'COMPLETE_WITH_ERRORS',tail(indices,1));cat('Worker finished.\n')
