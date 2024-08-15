## include generics


#'  Norm of a coordinate
#'
#' coord -> real
#' @export
Norm <- function(ob, ...) {
    UseMethod('Norm')
}


#'  Squared norm of a coordinate
#'
#' coord -> real
#' @export
Norm2 <- function(ob, ...) {
    UseMethod('Norm2')
}

#' Convert eligible objects to coords via generic S3 interface.
#' @export
as_coord <- function(ob, ...){
    UseMethod('as_coord')
}


#' Convert eligible objects to spatstat im type
#'
#' Will replace this eventually with an alternative implementation
#' @export
as.im <- function(...){
    UseMethod('as.im')
}


#' Prepare
#'
#' generic interface for post initialization preparation of spatial covariates
#' @export
prepare <- function(object,  dimyx, fractional){
  UseMethod('prepare')
}

#' @export
convert_raster <- function(object, ...){
  UseMethod('convert_raster')
}
