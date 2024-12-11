## Raster Type Covariates

#'  Rcov handles data passed in raster formats.
#' Currently, only 'im' type objects are handled
#' Rcov does minimal handling of the raster objects, its primary goal is for metadata management and allowing natural manipulation of raster valued covariates

#' @param interpolator
#' options include 'bin' which treats the raster as a piecewise constant surface given by the raster values,
#' and 'bilin' which implements a bilinear interpolation based on assigning the weights of the pixels to their upper left coordinate corners
#' alternatively, with the argument centroid = TRUE bilinear interpolation will be computed on the basis of centroids.
#'depending on how data is passed in, Rcov requires different types of arguments to correctly harmonize the data with other covariate representations
#' @export
Rcov <- function(..., dimyx = c(128, 128), W = NULL, interpolator = c('bilin', 'bin', 'spatstat', 'nearest')) {
    interpolator <- match.arg(interpolator, interpolator)
    rasts <- list2(...)

    prepped <- lapply(rasts, \(ob) Rcov_prepare(ob, dimyx = dimyx, W = W))

    vctrs::new_vctr(prepped,
                    interpolator = interpolator,
                    class = c('Rcov', 'spatial_covariate'))
}

## Rcov_impl <- function (ob, ...) {
##     stop()
## }

#' @exportS3Method
Rcov_prepare.terra <- function(){
    stop()
}

#' @exportS3Method
Rcov_prepare.raster <- function ()  {
    stop()
}


Rcov_prepare.im <- function (ob, ...) {
    ob
}

Rcov_prepare.stars <- function ( )
{
    stop()
}

#' @export
evaluate.Rcov <- function(object, locations, ...){
    stopifnot(is_coord(locations))

    x <- coordx(locations)
    y <- coordy(locations)

    interpolation_choice <- attr(object, 'interpolator')
    if (interpolation_choice == "bilin") {
        interpolator <- function(r) spatstat.geom::interp.im(r, x,y, bilinear = TRUE)

    }
     if (interpolation_choice == "spatstat") {   
         interpolator <-  function(r) spatstat.geom::interp.im(r, x,y, bilinear = FALSE)
     }

    do.call(c, lapply(object, interpolator))
}

#' @export
format.Rcov <- function(r, ...){
    rep('Raster Object', vec_size(r))
}

#' @export
vec_ptype_abbr.Rcov <- function(r, ...){
    'Rcov'
}


#' @export
Rcov_prepare <- function(object, ...) {
    UseMethod('Rcov_prepare')
}
#' @export
Rcov_prepare.im <- function(object, ...){
    return(object)
}

#' @export
Rcov_prepare.default <- function(object, ...){
    cl <- class(object)[[1]]
    message <- paste0('Objects of type ', cl, ' not currently supported.')
    stop(message)
}




if (FALSE){ 
rasterize_road <- function(road, nr = 512, nc = 512) {
    ovec <- vect(road)
    orast <- rast(ovec, nr, nc)

    linerast <- rasterize(ovec, orast, background = 0)
    matrast <- as.matrix(linerast, wide = TRUE)

    for (i in 1:ncol(matrast)){
        matrast[,i] <- rev(matrast[,i])
    }
    spatstat.geom::im(matrast, xrange = c(bb$xmin, bb$xmax), yrange = c(bb$ymin, bb$ymax))
}
}
