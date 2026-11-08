REPO_ROOT <- normalizePath(getwd(), winslash="/", mustWork=TRUE)
if (!file.exists(file.path(REPO_ROOT,"config","jsd.R"))) stop("Run from the repository root")
DATA_ROOT <- Sys.getenv("ADNI_DATA_DIR", file.path(REPO_ROOT,"data"))
RUN_ROOT <- Sys.getenv("ADNI_WORK_DIR", file.path(REPO_ROOT,"analysis"))
stage_dir <- function(name) {
  p <- file.path(RUN_ROOT,name)
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
  normalizePath(p,winslash="/",mustWork=TRUE)
}
source(file.path(REPO_ROOT,"src","json_io.R"))
