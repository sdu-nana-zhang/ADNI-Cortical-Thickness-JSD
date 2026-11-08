source('src/project.R')
suppressPackageStartupMessages(library(nlme))
root<-stage_dir('figures')
dir.create(file.path(root,'data'),showWarnings=FALSE)
a<-readRDS(file.path(stage_dir('core'),'results','fit_TRAJECTORY_ALL_AVAILABLE_L0.rds'))
d <- a$dat; pp <- d[!duplicated(d$PTID),]
r <- d[1,,drop=FALSE]
for(v in c('AGE_C','MALE_C','EDUC_C')) r[[v]] <- mean(pp[[v]])
out <- list()
for(g in c('CN','MCI','AD')) {
  rr <- r[rep(1,151),,drop=FALSE]
  rr$TIME <- seq(0,max(d$TIME[d$BASELINE_DX==g]),length.out=151)
  rr$BASELINE_DX <- factor(g,levels=c('CN','MCI','AD'))
  X <- model.matrix(delete.response(terms(a$form)),rr)[,names(a$b),drop=FALSE]
  est <- as.vector(X %*% a$b)
  se <- sqrt(rowSums((X %*% a$V)*X));q <- qt(.975,a$df)
  out[[g]] <- data.frame(GROUP=g,TIME=rr$TIME,ESTIMATE=est,CI_LOW=est-q*se,CI_HIGH=est+q*se)
}
write.csv(do.call(rbind,out),file.path(root,'data/adni_adjusted_curves.csv'),row.names=FALSE)
write.csv(data.frame(COVARIATE=c('AGE_C','MALE_C','EDUC_C'),VALUE=unlist(r[1,c('AGE_C','MALE_C','EDUC_C')])),file.path(root,'data/adni_prediction_reference.csv'),row.names=FALSE)
cat('Exported fixed-effect curves using stored coefficients and CR1 covariance; no model fitted.\n')
