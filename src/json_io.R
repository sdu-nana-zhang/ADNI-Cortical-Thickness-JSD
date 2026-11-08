fromJSON <- function(path) dget(sub("\\.json$", ".R", path))
json_encode <- function(x) {
 if(is.null(x)) return("null")
 if(is.list(x)) {
  if(!is.null(names(x))) return(paste0("{",paste(paste0(vapply(names(x),function(s)encodeString(s,quote='"'),character(1)),":",vapply(x,json_encode,character(1))),collapse=","),"}"))
  return(paste0("[",paste(vapply(x,json_encode,character(1)),collapse=","),"]"))
 }
 if(length(x)!=1L) return(paste0("[",paste(vapply(as.list(x),json_encode,character(1)),collapse=","),"]"))
 if(is.na(x)) return("null")
 if(is.logical(x)) return(if(x)"true" else "false")
 if(is.numeric(x)) return(if(is.finite(x))format(x,digits=17,scientific=FALSE,trim=TRUE) else "null")
 encodeString(as.character(x),quote='"')
}
toJSON <- function(x,...) json_encode(x)
write_json <- function(x,path,...) writeLines(json_encode(x),path,useBytes=TRUE)
