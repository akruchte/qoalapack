                                        # convenience functions for mgcv and pretty printing tools
## TODO
#' @exportS3Method
print.smooth.spec <- function(object, ...) {
    cat(paste0(object$label, " with basis class = ", class(object), ".\n"))
}

print.smooth <- function(object, ...) {
    cat(paste0("A smooth term.\n"))
}



#' @exportS3Method
print.tp.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.bs.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.ps.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.cr.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.te.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.ti.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.t2.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.t2.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.ds.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.cs.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.cc.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.sos.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.cp.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.re.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.mrf.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.gp.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.so.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.sw.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.sf.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.ad.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.sz.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}
#' @exportS3Method
print.fs.smooth.spec <- function(object, ... ) {print.smooth.spec(object, ...)}










## gratia package
