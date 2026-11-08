source("src/project.R")
suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8"))

suppressPackageStartupMessages(library(glmnet))

root <- stage_dir('score')
args <- commandArgs(trailingOnly=TRUE)
is_worker <- length(args)>0 && args[1]=="worker"
collect_only <- length(args)>0 && args[1]=="collect"
worker_id <- if(is_worker)as.integer(args[2]) else 1L
worker_count <- if(is_worker)as.integer(args[3]) else 1L
stopifnot(worker_id>=1,worker_id<=worker_count)
progress_name <- if(is_worker)sprintf("fit_worker%02d.log",worker_id) else "fit_progress.log"
state_name <- if(is_worker)sprintf("progress_worker%02d.json",worker_id) else "progress.json"
cfg <- fromJSON(file.path(root,"protocol","frozen_protocol.json"))
stopifnot(isTRUE(cfg$frozen_before_fit), cfg$n>1)
bounded_recovery <- is_worker && length(args)>=4 && args[4]=="bounded"
if(bounded_recovery){
 cfg$glmnet_maxit<-100000L
 cfg$glmnet_thresholds<-cfg$glmnet_thresholds[1]
 cat("Numerical recovery: bounded glmnet iterations followed by same-objective L-BFGS-B if needed.\n")
}
implementation_id <- "correct_population_sd_division_v2"
protocol_md5 <- unname(tools::md5sum(file.path(root,"protocol","frozen_protocol.json")))
input <- list.files(file.path(root,"inputs"),pattern="^classifier_input[.]csv$",full.names=TRUE)
stopifnot(length(input)==1)
d <- read.csv(input,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM")
stopifnot(nrow(d)==cfg$n,!anyDuplicated(d$PTID),all(d$Y %in% c(0,1)),all(table(d$Y)>0))
feature_names <- names(d)[6:ncol(d)]
X <- as.matrix(d[,feature_names]); storage.mode(X)<-"double"; y<-d$Y
stopifnot(ncol(X)==2020,all(is.finite(X)))
outer <- read.csv(file.path(root,"protocol","outer_patient_folds.csv"),fileEncoding="UTF-8-BOM")
inner <- read.csv(file.path(root,"protocol","inner_patient_folds.csv"),fileEncoding="UTF-8-BOM")
lt <- sort(cfg$traditional_penalties,decreasing=TRUE)
ratios <- sort(cfg$JSD_penalty_multipliers,decreasing=TRUE)
base_map <- list(M0=c("AGE","MALE","EDUCATION"),M2=c("AGE","MALE","EDUCATION","CORTEX_CM3","DKT62_AREA_WEIGHTED_THICKNESS_MM"))
base_map$M4 <- c(base_map$M2,grep("^ROI_",feature_names,value=TRUE))
jsd_names <- grep("^JSD_",feature_names,value=TRUE)
stopifnot(length(jsd_names)==1891,length(base_map$M4)==129)
pairs <- list(c("M0","M1"),c("M2","M3"),c("M4","M5"))
models <- paste0("M",0:5)
ng <- length(lt)
grid_base <- data.frame(LAMBDA_T=lt,RATIO=Inf,LAMBDA_J=Inf)
grid_jsd <- do.call(rbind,c(list(grid_base),lapply(ratios,function(r)data.frame(LAMBDA_T=lt,RATIO=r,LAMBDA_J=lt*r))))
grids <- setNames(lapply(models,function(m)if(m %in% c("M0","M2","M4"))grid_base else grid_jsd),models)
solver_thresholds <- unique(c(cfg$glmnet_thresholds,1e-16,1e-18))
glmnet.control(fdev=0,devmax=1)
logmsg <- function(s){cat(s,"\n");flush.console();cat(s,"\n",file=file.path(root,"work",progress_name),append=TRUE)}
prepare <- function(tr,va){
 a<-X[tr,,drop=FALSE];mu<-colMeans(a);sd<-sqrt(colMeans(sweep(a,2,mu,FUN="-")^2))
 keep<-is.finite(sd)&sd>.Machine$double.eps*pmax(1,abs(mu))
 stopifnot(all(is.finite(mu)),all(is.finite(sd)),sum(keep)>0)
 za<-sweep(sweep(a[,keep,drop=FALSE],2,mu[keep],FUN="-"),2,sd[keep],FUN="/")
 zv<-sweep(sweep(X[va,keep,drop=FALSE],2,mu[keep],FUN="-"),2,sd[keep],FUN="/")
 stopifnot(max(abs(colMeans(za)))<1e-10,max(abs(colMeans(za^2)-1))<1e-10)
 list(x=za,xt=zv,mean=mu,scale=sd,keep=keep)
}
fit_path <- function(z,ztest,ytr,traditional,ratio,tag){
 cols <- c(intersect(traditional,colnames(z)),if(is.finite(ratio))intersect(jsd_names,colnames(z)) else character())
 zi<-z[,cols,drop=FALSE];zv<-ztest[,cols,drop=FALSE]
 pen<-c(rep(1,sum(cols %in% traditional)),rep(ratio,sum(cols %in% jsd_names)))
 stopifnot(all(is.finite(pen)),all(pen>0),length(pen)==ncol(zi))
 fac<-sum(pen)/length(pen)
 design<-cbind(1,zi)
 checkgradient<-function(beta){
  pp<-plogis(design%*%beta)
  gr<-crossprod(design,sweep(pp,1,ytr))/nrow(zi)
  gr[-1,]<-gr[-1,,drop=FALSE]+beta[-1,,drop=FALSE]*(pen%o%lt)
  apply(abs(gr),2,max)
 }
 last<-NULL;good<-FALSE;warnings_all<-character()
 for(attempt in seq_along(cfg$glmnet_thresholds)){
  warns<-character()
  fit<-withCallingHandlers(glmnet(x=zi,y=ytr,family="binomial",alpha=0,lambda=lt*fac,penalty.factor=pen,standardize=FALSE,intercept=TRUE,thresh=cfg$glmnet_thresholds[attempt],maxit=cfg$glmnet_maxit),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart("muffleWarning")})
  warnings_all<-c(warnings_all,warns)
  if(length(fit$lambda)!=ng || fit$jerr!=0){
   warnings_all<-c(warnings_all,paste("glmnet incomplete path",length(fit$lambda),"jerr",fit$jerr,"threshold",cfg$glmnet_thresholds[attempt]))
   next
  }
  stopifnot(max(abs(fit$lambda-lt*fac))<1e-10)
  beta<-as.matrix(coef(fit));gin<-checkgradient(beta)
  if(all(gin<=cfg$optimizer_gradient_tolerance)){good<-TRUE;break}
 }
 used_fallback<-FALSE
 if(!good){
  used_fallback<-TRUE
  beta<-matrix(0,ncol(design),ng,dimnames=list(c("(Intercept)",cols),NULL))
  beta[1,]<-qlogis(mean(ytr));oldcoef<-as.matrix(coef(fit))
  for(j in seq_len(ng)){
   ix<-which(abs(fit$lambda-lt[j]*fac)<1e-10)
   if(length(ix))beta[,j]<-oldcoef[,ix[1]] else if(j>1)beta[,j]<-beta[,j-1]
   lp<-pen*lt[j]
   fn<-function(b){v<-as.numeric(design%*%b);mean(pmax(v,0)+log1p(exp(-abs(v)))-ytr*v)+.5*sum(lp*b[-1]^2)}
   gr<-function(b){v<-as.numeric(crossprod(design,plogis(design%*%b)-ytr))/nrow(design);v[-1]<-v[-1]+lp*b[-1];v}
   g0<-max(abs(gr(beta[,j])))
   if(g0>cfg$optimizer_gradient_tolerance){
    opt<-optim(beta[,j],fn,gr,method="L-BFGS-B",control=list(maxit=20000,factr=0,pgtol=1e-8,lmm=15))
    beta[,j]<-opt$par
    if(max(abs(gr(beta[,j])))>cfg$optimizer_gradient_tolerance){
     opt<-optim(beta[,j],fn,gr,method="L-BFGS-B",control=list(maxit=40000,factr=0,pgtol=1e-10,lmm=15))
     beta[,j]<-opt$par
    }
    warnings_all<-c(warnings_all,paste("L-BFGS-B numerical polish lambdaT",lt[j],"gradient",max(abs(gr(beta[,j]))),"convergence",opt$convergence,"function_evaluations",opt$counts[1]))
   }
  }
  gin<-checkgradient(beta)
  logmsg(paste("Numerical L-BFGS-B fallback",tag,"max gradient",max(gin)))
 }
 stopifnot(all(gin<=cfg$optimizer_gradient_tolerance),all(is.finite(beta)))
 pred<-plogis(cbind(1,zv)%*%beta)
 if(!used_fallback){reference<-as.matrix(predict(fit,newx=zv,type="response"));stopifnot(max(abs(pred-reference))<1e-10)}
 stopifnot(all(is.finite(pred)))
 di<-data.frame(TAG=tag,LAMBDA_T=lt,LAMBDA_J=if(is.finite(ratio))lt*ratio else Inf,RATIO=ratio,P=ncol(zi),GRAD_INF=gin,CONVERGED=TRUE,ATTEMPTS=attempt+as.integer(used_fallback),NPASSES=fit$npasses,JERR=0,WARNING=if(length(warnings_all))paste(unique(warnings_all),collapse="; ") else "")
 list(pred=pred,beta=beta,cols=cols,pen=pen,diagnostics=di)
}
choose <- function(loss,grid){
 ix<-which(loss<=min(loss)+1e-10*max(1,abs(min(loss))))
 q<-grid[ix,,drop=FALSE]
 ix[order(!is.infinite(q$LAMBDA_J),-q$LAMBDA_J,-q$LAMBDA_T)[1]]
}
all_results<-list();started<-proc.time()[3]
for(re in 1:5){
 of<-outer[outer$REPEAT==re,];stopifnot(identical(as.character(of$PTID),as.character(d$PTID)),identical(as.integer(of$ROW1),seq_len(nrow(d))))
 for(fo in 1:5){
  if(is_worker && (((re-1)*5+fo-1) %% worker_count + 1L)!=worker_id)next
  key<-sprintf("repeat%02d_fold%02d",re,fo);checkpoint<-file.path(root,"models",paste0(key,".rds"))
  if(file.exists(checkpoint)){
   obj<-readRDS(checkpoint);stopifnot(identical(obj$protocol_md5,protocol_md5),identical(obj$implementation_id,implementation_id));all_results[[key]]<-obj;logmsg(paste("Reused verified checkpoint",key));next
  }
  if(collect_only)stop(paste("Missing checkpoint",key))
  tr<-which(of$OUTER_FOLD!=fo);te<-which(of$OUTER_FOLD==fo)
  stopifnot(length(intersect(d$PTID[tr],d$PTID[te]))==0)
  inf<-inner[inner$REPEAT==re & inner$OUTER_FOLD==fo,];stopifnot(identical(as.integer(inf$ROW1),as.integer(tr)),identical(as.character(inf$PTID),as.character(d$PTID[tr])))
  ip<-setNames(lapply(models,function(m)matrix(NA_real_,length(tr),nrow(grids[[m]]))),models)
  states<-list();diagnostics<-list();audit_rows<-list();nd<-0
  for(fi in 1:4){
   it<-tr[inf$INNER_FOLD!=fi];iv<-tr[inf$INNER_FOLD==fi];local<-which(inf$INNER_FOLD==fi)
   stopifnot(!length(intersect(d$PTID[it],d$PTID[iv])),!length(intersect(d$PTID[it],d$PTID[te])),!length(intersect(d$PTID[iv],d$PTID[te])))
   st<-prepare(it,iv);states[[as.character(fi)]]<-list(train=it,validation=iv,mean=st$mean,scale=st$scale,keep=st$keep)
   for(pair in pairs){
    bm<-pair[1];jm<-pair[2]
    fitb<-fit_path(st$x,st$xt,y[it],base_map[[bm]],Inf,paste(key,"inner",fi,bm))
    ip[[bm]][local,]<-fitb$pred;ip[[jm]][local,1:ng]<-fitb$pred
    nd<-nd+1;diagnostics[[nd]]<-fitb$diagnostics
    for(k in seq_along(ratios)){
     fitj<-fit_path(st$x,st$xt,y[it],base_map[[bm]],ratios[k],paste(key,"inner",fi,jm,"ratio",ratios[k]))
     ip[[jm]][local,ng+(k-1)*ng+1:ng]<-fitj$pred
     nd<-nd+1;diagnostics[[nd]]<-fitj$diagnostics
    }
   }
   audit_rows[[fi]]<-data.frame(REPEAT=re,OUTER_FOLD=fo,INNER_FOLD=fi,N_TRAIN=length(it),N_VALIDATION=length(iv),N_OUTER_TEST=length(te),N_OVERLAP=0)
   logmsg(sprintf("%s inner %d/4 complete; elapsed %.1f s",key,fi,proc.time()[3]-started))
  }
  losses<-list();selection<-list();allgrid<-list()
  for(m in models){
   stopifnot(all(is.finite(ip[[m]])))
   p<-pmin(pmax(ip[[m]],1e-15),1-1e-15)
   loss<-colMeans(-y[tr]*log(p)-(1-y[tr])*log1p(-p))
   losses[[m]]<-loss;sel<-choose(loss,grids[[m]]);selection[[m]]<-sel
   gg<-grids[[m]];gg$MODEL<-m;gg$REPEAT<-re;gg$OUTER_FOLD<-fo;gg$CANDIDATE<-seq_len(nrow(gg));gg$INNER_LOGLOSS<-loss;gg$SELECTED<-as.integer(seq_len(nrow(gg))==sel);allgrid[[m]]<-gg
  }
  for(pair in pairs)stopifnot(identical(ip[[pair[1]]],ip[[pair[2]]][,1:ng,drop=FALSE]))
  st<-prepare(tr,te);pred<-matrix(NA_real_,length(te),6,dimnames=list(NULL,models));frozen<-list();params<-list()
  for(pair in pairs){
   bm<-pair[1];jm<-pair[2]
   fitb<-fit_path(st$x,st$xt,y[tr],base_map[[bm]],Inf,paste(key,"outer",bm));nd<-nd+1;diagnostics[[nd]]<-fitb$diagnostics
   for(m in pair){
    sel<-selection[[m]];g<-grids[[m]][sel,]
    if(is.infinite(g$RATIO))fit<-fitb else{
     fit<-fit_path(st$x,st$xt,y[tr],base_map[[bm]],g$RATIO,paste(key,"outer",m));nd<-nd+1;diagnostics[[nd]]<-fit$diagnostics
    }
    li<-which(lt==g$LAMBDA_T);stopifnot(length(li)==1)
    pred[,m]<-fit$pred[,li]
    original_index<-match(fit$cols,feature_names)
    frozen[[m]]<-list(beta=fit$beta[,li,drop=FALSE],features=fit$cols,original_indices=original_index,mean=st$mean[original_index],scale=st$scale[original_index],penalty=fit$pen*g$LAMBDA_T,lambdaT=g$LAMBDA_T,lambdaJ=g$LAMBDA_J,training_PTID=d$PTID[tr],test_PTID=d$PTID[te])
    params[[m]]<-data.frame(REPEAT=re,OUTER_FOLD=fo,MODEL=m,TRAIN_N=length(tr),TEST_N=length(te),LAMBDA_T=g$LAMBDA_T,LAMBDA_J=g$LAMBDA_J,RATIO=g$RATIO,INNER_LOGLOSS=losses[[m]][sel],P=length(fit$cols),GRAD_INF=fit$diagnostics$GRAD_INF[li],JSD_OFF=is.infinite(g$RATIO))
   }
  }
  stopifnot(all(is.finite(pred)),all(pred>=0),all(pred<=1))
  obj<-list(implementation_id=implementation_id,protocol_md5=protocol_md5,repeat_id=re,outer_fold=fo,train=tr,test=te,outer_preprocessing=list(mean=st$mean,scale=st$scale,keep=st$keep),inner_preprocessing=states,inner_predictions=ip,inner_grid=do.call(rbind,allgrid),parameters=do.call(rbind,params),diagnostics=do.call(rbind,diagnostics),partition_audit=do.call(rbind,audit_rows),models=frozen,predictions=pred)
  saveRDS(obj,checkpoint,compress="gzip");all_results[[key]]<-obj
  logmsg(sprintf("%s OUTER COMPLETE; elapsed %.1f s",key,proc.time()[3]-started))
  write_json(list(status="RUNNING",completed_outer=length(all_results),total_outer=25,elapsed_seconds=proc.time()[3]-started),file.path(root,"work",state_name),auto_unbox=TRUE,pretty=TRUE)
 }
}
if(is_worker){
 write_json(list(status="COMPLETE",worker=worker_id,completed_outer=length(all_results)),file.path(root,"work",state_name),auto_unbox=TRUE,pretty=TRUE)
 logmsg(paste("WORKER COMPLETE",worker_id));quit(save="no",status=0)
}
stopifnot(length(all_results)==25)
rows<-list();k<-0
for(o in all_results)for(m in models){k<-k+1;rows[[k]]<-data.frame(PTID=d$PTID[o$test],ROW1=o$test,CLASS=d$CLASS[o$test],Y=y[o$test],REPEAT=o$repeat_id,OUTER_FOLD=o$outer_fold,MODEL=m,P_AD=o$predictions[,m])}
write.csv(do.call(rbind,rows),file.path(root,"results","OOF_predictions_long.csv"),row.names=FALSE,fileEncoding="UTF-8")
write.csv(do.call(rbind,lapply(all_results,`[[`,"parameters")),file.path(root,"results","outer_selected_parameters.csv"),row.names=FALSE)
write.csv(do.call(rbind,lapply(all_results,`[[`,"inner_grid")),file.path(root,"results","inner_grid_losses.csv"),row.names=FALSE)
write.csv(do.call(rbind,lapply(all_results,`[[`,"diagnostics")),file.path(root,"audit","optimizer_diagnostics.csv"),row.names=FALSE)
write.csv(do.call(rbind,lapply(all_results,`[[`,"partition_audit")),file.path(root,"audit","actual_partition_overlap.csv"),row.names=FALSE)
for(m in models){
 o<-all_results[[1]];mm<-o$models[[m]]
 tab<-data.frame(FEATURE=c("INTERCEPT",mm$features),BETA=as.numeric(mm$beta),MEAN=c(0,mm$mean),SCALE=c(1,mm$scale),PENALTY=c(0,mm$penalty))
 write.csv(tab,file.path(root,"audit",paste0("numerical_check_",m,".csv")),row.names=FALSE)
}
all_diag<-do.call(rbind,lapply(all_results,`[[`,"diagnostics"))
record<-list(status="COMPLETE",n=nrow(d),class_counts=list(CN=sum(y==0),AD=sum(y==1)),outer_repeats=5,outer_folds=5,inner_folds=4,n_outer_predictions=nrow(do.call(rbind,rows)),n_optimizer_candidate_checks=nrow(all_diag),max_gradient=max(all_diag$GRAD_INF),all_optimizer_converged=all(all_diag$CONVERGED),n_numerical_retry_paths=length(unique(all_diag$TAG[all_diag$ATTEMPTS>1])),protocol_md5=protocol_md5,implementation_id=implementation_id, solver_thresholds=1e-12, glmnet_iteration_budget=100000L, numerical_fallback="stats::optim L-BFGS-B; same objective, grid and gradient tolerance", parallel_worker_count=if(collect_only)NA_integer_ else worker_count, glmnet_version=as.character(packageVersion("glmnet")),R_version=R.version.string,collection_elapsed_seconds=proc.time()[3]-started,imputation=FALSE,no_future_outcome_features=TRUE)
write_json(record,file.path(root,"audit","execution_record.json"),auto_unbox=TRUE,pretty=TRUE)
write_json(list(status="COMPLETE",completed_outer=25,total_outer=25),file.path(root,"work",state_name),auto_unbox=TRUE,pretty=TRUE)
capture.output(sessionInfo(),file=file.path(root,"audit","R_sessionInfo.txt"))
logmsg(toJSON(record,auto_unbox=TRUE))
