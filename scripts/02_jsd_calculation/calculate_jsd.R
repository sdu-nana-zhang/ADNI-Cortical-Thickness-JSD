source("src/project.R")
source("src/surface.R")
source("src/jsd.R")
cfg <- thickness_config()
manifest <- read.csv(file.path(DATA_ROOT,"scan_manifest.csv"),check.names=FALSE,stringsAsFactors=FALSE)
out <- stage_dir("jsd")
dir.create(file.path(out,"scans"),showWarnings=FALSE)
rows <- list(); qc <- list(); stats <- list()
for(i in seq_len(nrow(manifest))) {
  a <- readRDS(file.path(stage_dir("roi"),sprintf("scan_%06d.rds",i)))
  z <- thickness_features(a$raw,cfg,return_mass=TRUE)
  if(nrow(z$roi_stats)!=62L || !all(z$qc$TechnicalPass) || any(!is.finite(z$edges$JSD)))
    stop("ROI or JSD QC failed at manifest row ",i,"; review the authorized input manifest")
  mu <- z$roi_stats$Mean[match(cfg$roi_names,z$roi_stats$ROI)]
  areas <- a$cortical$roi$AREA_MM2[match(cfg$roi_names,a$cortical$roi$ROI)]
  avg <- if(all(is.finite(c(mu,areas)))&&sum(areas)>0)
    global_mean_thickness(matrix(mu,nrow=1),matrix(areas,nrow=1)) else NA_real_
  values <- c(setNames(mu,paste0("thickness__Mean__",cfg$roi_names)),
    setNames(z$roi_stats$SD[match(cfg$roi_names,z$roi_stats$ROI)],paste0("thickness__SD__",cfg$roi_names)),
    setNames(z$edges$JSD,z$edges$FEATURE),CORTEX_CM3=a$cortical$cortex_volume_mm3/1000,
    DKT62_AREA_WEIGHTED_THICKNESS_MM=avg)
  rows[[i]] <- cbind(a$metadata,as.data.frame(as.list(values),check.names=FALSE))
  z$metadata <- a$metadata
  saveRDS(z,file.path(out,"scans",sprintf("scan_%06d.rds",i)),compress="gzip")
  qc[[i]] <- cbind(ROW=i,z$qc)
  stats[[i]] <- cbind(ROW=i,a$cortical$roi)
  message("JSD scan ",i," / ",nrow(manifest))
}
write.csv(do.call(rbind,rows),file.path(out,"JSD_features.csv"),row.names=FALSE,fileEncoding="UTF-8")
write.csv(do.call(rbind,qc),file.path(out,"ROI_QC.csv"),row.names=FALSE,fileEncoding="UTF-8")
write.csv(do.call(rbind,stats),file.path(out,"cortical_stats.csv"),row.names=FALSE,fileEncoding="UTF-8")
write.csv(edge_dictionary(),file.path(out,"edge_order.csv"),row.names=FALSE)
