#' @importFrom spatstat.geom as.im
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


#' @export
convert_raster <- function(object, ...){
  UseMethod('convert_raster')
}
