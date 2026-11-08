options(stringsAsFactors=FALSE)
source("src/project.R")
suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8"))
suppressPackageStartupMessages(library(nlme))
root <- stage_dir("core")
for(nm in c("results","audit")) dir.create(file.path(root,nm),showWarnings=FALSE)
read_csv <- function(path) read.csv(path, fileEncoding="UTF-8-BOM", check.names=FALSE, stringsAsFactors=FALSE)
write_csv <- function(x,path) write.csv(x,path,row.names=FALSE,fileEncoding="UTF-8",na="")
d <- read_csv(file.path(stage_dir("score"),"results","scored_MRI.csv"))
d <- d[,c("PTID","MRI_IMAGE_ID","MRI_DATE","DX_STD","AGE","MALE","EDUCATION","EXAMDATE","JSD_SCORE_Z","JSD_BASELINE_ROW","CORTEX_CM3","DKT62_AREA_WEIGHTED_THICKNESS_MM","ADAS13","MMSE")]
d <- d[order(d$PTID,d$MRI_DATE,d$MRI_IMAGE_ID),]
d$PTID <- as.character(d$PTID)
first <- d[!duplicated(d$PTID), c("PTID","MRI_DATE","DX_STD","AGE","MALE","EDUCATION")]
names(first) <- c("PTID","BASELINE_MRI_DATE","BASELINE_DX","BASELINE_AGE","BASELINE_MALE","BASELINE_EDUCATION")
d <- merge(d, first, by="PTID", all.x=TRUE, sort=FALSE)
d <- d[order(match(d$PTID,first$PTID),d$MRI_DATE,d$MRI_IMAGE_ID),]
d$BASELINE_SEX <- ifelse(d$BASELINE_MALE==1,"Male","Female")
d$BASELINE_DX <- factor(d$BASELINE_DX, levels=c("CN","MCI","AD"))
d$BASELINE_AGE <- as.numeric(d$BASELINE_AGE); d$BASELINE_EDUCATION <- as.numeric(d$BASELINE_EDUCATION)
d$YEARS_FROM_FIRST_MRI <- as.numeric(as.Date(d$MRI_DATE)-as.Date(d$BASELINE_MRI_DATE))/365.25
d$CLINICAL_YEARS_FROM_FIRST_MRI <- as.numeric(as.Date(d$EXAMDATE)-as.Date(d$BASELINE_MRI_DATE))/365.25
d$AGE_C <- d$BASELINE_AGE-75; d$MALE_C <- d$BASELINE_MALE; d$EDUC_C <- d$BASELINE_EDUCATION-16
d$CORTEX100 <- d$CORTEX_CM3/100; d$THICKNESS <- d$DKT62_AREA_WEIGHTED_THICKNESS_MM
fit_model <- function(dat,form,label) {
 chosen<-NULL;typ<-NULL
 for(kind in c("DIAGONAL_INTERCEPT_SLOPE","INTERCEPT_ONLY")) { for(opt in c("nlminb","optim")) {
  ww<-character(); ans<-tryCatch(withCallingHandlers(lme(fixed=form,data=dat,random=if(kind=="DIAGONAL_INTERCEPT_SLOPE") list(PTID=pdDiag(~TIME)) else ~1|PTID,method="REML",na.action=na.fail,control=lmeControl(opt=opt,maxIter=200,msMaxIter=200,niterEM=30,tolerance=1e-6,returnObject=FALSE)),warning=function(w){ww<<-c(ww,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e) NULL)
  if(!is.null(ans)&&all(is.finite(fixef(ans)))&&is.finite(ans$sigma)&&!any(grepl("converg|iteration limit",ww,ignore.case=TRUE))){chosen<-ans;typ<-kind;break}
 }; if(!is.null(chosen)) break }
 if(is.null(chosen)) stop(paste("No valid model",label))
 fit<-chosen;b<-fixef(fit);X<-model.matrix(delete.response(terms(form)),dat);Y<-model.response(model.frame(form,dat));D<-as.matrix(getVarCov(fit,type="random.effects"));sigma2<-fit$sigma^2;groups<-split(seq_len(nrow(dat)),dat$PTID,drop=TRUE);A<-matrix(0,ncol(X),ncol(X));U<-matrix(0,length(groups),ncol(X));rhs<-numeric(ncol(X));uix<-0
 for(ix in groups){Z<-if(typ=="DIAGONAL_INTERCEPT_SLOPE") cbind(1,dat$TIME[ix]) else matrix(1,length(ix),1);V<-Z%*%D%*%t(Z)+diag(sigma2,length(ix));W<-solve(V);xx<-X[ix,,drop=FALSE];A<-A+crossprod(xx,W%*%xx);rhs<-rhs+as.numeric(crossprod(xx,W%*%Y[ix]));uix<-uix+1;U[uix,]<-as.numeric(crossprod(xx,W%*%(Y[ix]-as.vector(xx%*%b))));}
 bread<-solve(A);G<-length(groups);N<-nrow(dat);P<-ncol(X);cr1<-G/(G-1)*(N-1)/(N-P);Vrob<-cr1*bread%*%crossprod(U)%*%bread;dimnames(Vrob)<-list(names(b),names(b));se<-sqrt(diag(Vrob));df<-G-1;crit<-qt(.975,df)
 co<-data.frame(ANALYSIS=label,TERM=names(b),BETA=as.numeric(b),CR1_SE=se,CI_LOW=b-crit*se,CI_HIGH=b+crit*se,P=2*pt(-abs(b/se),df),DF=df,N_PATIENTS=G,N_MRI=N,MODEL_BASED_SE=sqrt(diag(fit$varFix)))
 list(fit=fit,b=b,V=Vrob,df=df,dat=dat,form=form,label=label,coef=co,random_type=typ)
}
contrast <- function(ans,L,type,name){L<-as.numeric(L);est<-sum(L*ans$b);se<-sqrt(as.numeric(t(L)%*%ans$V%*%L));crit<-qt(.975,ans$df);data.frame(ANALYSIS=ans$label,TYPE=type,CONTRAST=name,ESTIMATE=est,CR1_SE=se,CI_LOW=est-crit*se,CI_HIGH=est+crit*se,P=2*pt(-abs(est/se),ans$df),DF=ans$df,N_PATIENTS=length(unique(ans$dat$PTID)),N_MRI=nrow(ans$dat))}
reference <- function(ans){x<-ans$dat;pp<-x[!duplicated(x$PTID),];r<-x[1,,drop=FALSE];r$AGE_C<-mean(pp$AGE_C);r$MALE_C<-mean(pp$MALE_C);r$EDUC_C<-mean(pp$EDUC_C);r$CORTEX_BETWEEN<-mean(pp$CORTEX_BETWEEN);r$THICKNESS_BETWEEN<-mean(pp$THICKNESS_BETWEEN);r$CORTEX_WITHIN<-0;r$THICKNESS_WITHIN<-0;r}
mm <- function(ans,ref,time,dx){ref$TIME<-time;ref$BASELINE_DX<-factor(dx,levels=c("CN","MCI","AD"));model.matrix(delete.response(terms(ans$form)),ref)[,names(ans$b),drop=FALSE]}
add_decomp <- function(x){for(v in c("JSD_SCORE_Z","CORTEX100","THICKNESS")){p<-c(JSD_SCORE_Z="SCORE",CORTEX100="CORTEX",THICKNESS="THICKNESS")[[v]];x[[paste0(p,"_BETWEEN")]]<-ave(x[[v]],x$PTID,FUN=mean);x[[paste0(p,"_WITHIN")]]<-x[[v]]-x[[paste0(p,"_BETWEEN")]]};x}
trajectory <- function(dat,label,with_mri){dat$TIME<-dat$YEARS_FROM_FIRST_MRI;dat<-add_decomp(dat);form<-JSD_SCORE_Z~TIME*BASELINE_DX+(AGE_C+MALE_C+EDUC_C)*TIME;if(with_mri) form<-update(form,.~.+CORTEX_BETWEEN+CORTEX_WITHIN+THICKNESS_BETWEEN+THICKNESS_WITHIN);a<-fit_model(dat,form,label);ref<-reference(a);slopes<-list();out<-list();for(g in c("CN","MCI","AD")){L<-mm(a,ref,1,g)-mm(a,ref,0,g);slopes[[g]]<-L;out[[length(out)+1]]<-contrast(a,L,"GROUP_SLOPE",g)};for(pair in list(c("MCI","CN"),c("AD","CN"),c("AD","MCI")))out[[length(out)+1]]<-contrast(a,slopes[[pair[1]]]-slopes[[pair[2]]],"SLOPE_DIFFERENCE",paste(pair,collapse="-"));list(model=a,results=do.call(rbind,out),coef=a$coef,sample=data.frame(ANALYSIS=label,BASELINE_DX=as.character(dat$BASELINE_DX),PTID=as.character(dat$PTID),MRI_IMAGE_ID=dat$MRI_IMAGE_ID,TIME=dat$TIME,JSD_SCORE_Z=dat$JSD_SCORE_Z))}
res<-list();coefs<-list();samps<-list();fits<-list()
core_req<-c("JSD_SCORE_Z","YEARS_FROM_FIRST_MRI","AGE_C","MALE_C","EDUC_C","CORTEX100","THICKNESS");print(colSums(is.na(d[,core_req])))
x<-d[complete.cases(d[,core_req]),];x<-x[as.numeric(ave(x$MRI_DATE,x$PTID,FUN=function(z) length(unique(z))))>=2,];x$PTID<-factor(x$PTID);cat("CORE_INPUT",nrow(x),length(unique(x$PTID)),"\n")
if(MODEL_SCOPE=="trajectory") for(w in c(FALSE,TRUE)){label<-if(w) "TRAJECTORY_ALL_AVAILABLE_L1" else "TRAJECTORY_ALL_AVAILABLE_L0";a<-trajectory(x,label,w);a$results$MODEL<-if(w) "L1" else "L0";a$results$WINDOW<-"ALL_AVAILABLE";res[[label]]<-a$results;coefs[[label]]<-a$coef;samps[[label]]<-a$sample;fits[[label]]<-a$model}
if(MODEL_SCOPE=="cognition") for(ep in c("ADAS13","MMSE")){req<-c(ep,"JSD_SCORE_Z","CLINICAL_YEARS_FROM_FIRST_MRI","AGE_C","MALE_C","EDUC_C","CORTEX100","THICKNESS");z<-d[complete.cases(d[,req]),];z<-z[as.numeric(ave(z$MRI_DATE,z$PTID,FUN=function(v) length(unique(v))))>=2,];z<-z[as.numeric(ave(z$EXAMDATE,z$PTID,FUN=function(v) length(unique(v))))>=2,];z$TIME<-z$CLINICAL_YEARS_FROM_FIRST_MRI;z$PTID<-factor(z$PTID);z$COGNITION<-z[[ep]];z<-add_decomp(z);form<-COGNITION~TIME*BASELINE_DX+(AGE_C+MALE_C+EDUC_C)*TIME+SCORE_BETWEEN+SCORE_WITHIN+CORTEX_BETWEEN+CORTEX_WITHIN+THICKNESS_BETWEEN+THICKNESS_WITHIN;label<-paste("COGNITION",ep,"ALL_AVAILABLE_C1",sep="_");a<-fit_model(z,form,label);q<-a$coef[a$coef$TERM%in%c("SCORE_BETWEEN","SCORE_WITHIN"),];q$OUTCOME<-ep;q$WINDOW<-"ALL_AVAILABLE";q$MODEL<-"C1";res[[label]]<-q;coefs[[label]]<-a$coef;samps[[label]]<-data.frame(ANALYSIS=label,PTID=as.character(z$PTID),MRI_IMAGE_ID=z$MRI_IMAGE_ID,TIME=z$TIME,COGNITION=z$COGNITION,SCORE_BETWEEN=z$SCORE_BETWEEN,SCORE_WITHIN=z$SCORE_WITHIN);fits[[label]]<-a$model}

write_csv(do.call(rbind,res),file.path(root,"results",paste0(MODEL_SCOPE,"_results.csv")))
write_csv(do.call(rbind,coefs),file.path(root,"results",paste0(MODEL_SCOPE,"_coefficients.csv")))
write_csv(do.call(rbind,samps),file.path(root,"audit",paste0(MODEL_SCOPE,"_included.csv")))
for(nm in names(fits)) saveRDS(fits[[nm]],file.path(root,"results",paste0("fit_",nm,".rds")),compress="gzip")
print(do.call(rbind,res),row.names=FALSE)
