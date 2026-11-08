source("src/project.R")
d<-read.csv(file.path(stage_dir("score"),"results","scored_MRI.csv"),fileEncoding="UTF-8-BOM",check.names=FALSE)
baseline<-read.csv(file.path(stage_dir("score"),"results","baseline_scores.csv"),fileEncoding="UTF-8-BOM",check.names=FALSE)
d<-merge(baseline[,c("PTID","MRI_IMAGE_ID")],d,by=c("PTID","MRI_IMAGE_ID"),sort=FALSE)
d$DX_STD<-factor(d$DX_STD,levels=c('CN','MCI','AD'))
d$AGE<-as.numeric(d$AGE);d$MALE<-as.numeric(d$MALE);d$EDUCATION<-as.numeric(d$EDUCATION)
d$CORTEX_CM3<-as.numeric(d$CORTEX_CM3);d$DKT62_AREA_WEIGHTED_THICKNESS_MM<-as.numeric(d$DKT62_AREA_WEIGHTED_THICKNESS_MM)
d$JSD_SCORE_Z<-as.numeric(d$JSD_SCORE_Z)
d$ADAS13_ANALYSIS<-as.numeric(d$ADAS13);d$MMSE_ANALYSIS<-as.numeric(d$MMSE)
hc3<-function(f,dat){dd<-dat[complete.cases(dat[,all.vars(f)]),,drop=FALSE];fit<-lm(f,data=dd);X<-model.matrix(fit);e<-resid(fit);h<-hatvalues(fit);bread<-solve(crossprod(X));V<-bread%*%crossprod(X*as.numeric(e/(1-h)))%*%bread;b<-coef(fit);se<-sqrt(diag(V));df<-df.residual(fit);q<-qt(.975,df);co<-data.frame(TERM=names(b),BETA=as.numeric(b),HC3_SE=se,CI_LOW=b-q*se,CI_HIGH=b+q*se,P=2*pt(-abs(b/se),df),DF=df,N=nrow(dd),CN=sum(dd$DX_STD=='CN'),MCI=sum(dd$DX_STD=='MCI'),AD=sum(dd$DX_STD=='AD'));list(co=co,n=nrow(dd))}
forms<-list(C0=~JSD_SCORE_Z+AGE+MALE+EDUCATION,C1=~JSD_SCORE_Z+AGE+MALE+EDUCATION+CORTEX_CM3+DKT62_AREA_WEIGHTED_THICKNESS_MM,C2=~JSD_SCORE_Z+AGE+MALE+EDUCATION+CORTEX_CM3+DKT62_AREA_WEIGHTED_THICKNESS_MM+DX_STD)
out<-list();k<-0
for(y in c('ADAS13_ANALYSIS','MMSE_ANALYSIS'))for(m in names(forms)){k<-k+1;f<-as.formula(paste(y,paste(deparse(forms[[m]]),collapse=' ')));a<-hc3(f,d);a$co$OUTCOME<-sub('_ANALYSIS','',y);a$co$MODEL<-m;out[[k]]<-a$co}
write.csv(do.call(rbind,out),file.path(stage_dir('core'),'results','baseline_coefficients.csv'),row.names=FALSE)
print(subset(do.call(rbind,out),TERM=='JSD_SCORE_Z'),row.names=FALSE)
mci_out<-list();z<-d[d$DX_STD=='MCI',,drop=FALSE];for(y in c('ADAS13_ANALYSIS','MMSE_ANALYSIS')){f<-as.formula(paste(y,paste(deparse(forms$C1),collapse=' ')));a<-hc3(f,z);a$co$OUTCOME<-sub('_ANALYSIS','',y);a$co$MODEL<-'C1_MCI_ONLY';mci_out[[y]]<-a$co;print(subset(a$co,TERM=='JSD_SCORE_Z'),row.names=FALSE)}

write.csv(do.call(rbind,mci_out),file.path(stage_dir("core"),"results","baseline_MCI_coefficients.csv"),row.names=FALSE)
