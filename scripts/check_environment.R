invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
required<-c(glmnet='5.1',nlme='3.1.169',Matrix='1.7.5')
for(p in names(required)) {
  if(!requireNamespace(p,quietly=TRUE))stop('Missing R package: ',p)
  v<-as.character(packageVersion(p))
  cat(p,v,if(v==required[[p]])'recorded version' else 'version differs from the validated runtime','\n')
}
cat(R.version.string,'\n')
