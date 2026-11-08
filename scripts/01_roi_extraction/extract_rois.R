source("src/project.R")
source("src/surface.R")
source("src/jsd.R")
args <- commandArgs(trailingOnly=TRUE)
worker <- if(length(args)) as.integer(args[1]) else 1L
workers <- if(length(args)>1) as.integer(args[2]) else 1L
stopifnot(worker>=1L, workers>=worker)
manifest <- read.csv(file.path(DATA_ROOT,"scan_manifest.csv"),check.names=FALSE,stringsAsFactors=FALSE)
required <- c("PTID","MRI_IMAGE_ID","MRI_DATE","FOLDER")
stopifnot(all(required %in% names(manifest)), !anyDuplicated(paste(manifest$PTID,manifest$MRI_IMAGE_ID)))
fsroot <- Sys.getenv("ADNI_FS_ROOT")
if(!nzchar(fsroot)) stop("Set ADNI_FS_ROOT to the reconstructed-subjects directory")
fsroot <- normalizePath(fsroot,winslash="/",mustWork=TRUE)
out <- stage_dir("roi")
indices <- seq(worker,nrow(manifest),by=workers)
for(i in indices) {
  path <- normalizePath(file.path(fsroot,manifest$FOLDER[i]),winslash="/",mustWork=TRUE)
  if(!startsWith(tolower(path),paste0(tolower(fsroot),"/"))) stop("Scan path escapes ADNI_FS_ROOT")
  raw <- surface_extract_scan(path)
  cortical <- cortical_mri_stats(path)
  saveRDS(list(raw=raw,cortical=cortical,metadata=manifest[i,,drop=FALSE]),
          file.path(out,sprintf("scan_%06d.rds",i)),compress="gzip")
  message("ROI scan ",i," / ",nrow(manifest))
}
