source('src/project.R')
root <- stage_dir('score')
files <- list.files(file.path(root,'models'),pattern='^repeat[0-9]{2}_fold[0-9]{2}\\.rds$',full.names=TRUE)
stopifnot(length(files)==25)
out <- list(); k <- 0L
for (f in sort(files)) {
 o <- readRDS(f); m <- o$models$M3; p <- o$parameters[o$parameters$MODEL=='M3',,drop=FALSE]
 rows <- data.frame(REPEAT=o$repeat_id,OUTER_FOLD=o$outer_fold,MODEL='M3',FEATURE=c('(Intercept)',m$features),BETA=as.numeric(m$beta),MEAN=c(0,m$mean),SCALE=c(1,m$scale),LAMBDA_T=p$LAMBDA_T,LAMBDA_J=p$LAMBDA_J,RATIO=p$RATIO,TRAIN_N=p$TRAIN_N,TEST_N=p$TEST_N,INNER_LOGLOSS=p$INNER_LOGLOSS,stringsAsFactors=FALSE)
 k<-k+1L; out[[k]]<-rows
}
par <- do.call(rbind,out)
write.csv(par,file.path(root,'results','M3_parameters.csv'),row.names=FALSE,fileEncoding='UTF-8')
cat('M3_PARAMETER_ROWS=',nrow(par),'\\n',sep='')
