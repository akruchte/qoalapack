## library(rlang)
## library(vctrs)
## library(spatstat.core)
## source('coord.R')
## source('covariate_types.R')

## Rcov currently uses im objects as the internal raster representation

## Rcov handles data passed in raster formats.
## initially, only 'im' type objects are handled
## Rcov does minimal handling of the raster objects,
## each raster object is assigned a personal interpolator function that determines how interpolation
## will be handled

## interpolator options include 'bin' which treats the raster as a piecewise constant surface given by the raster values,
## and 'bilin' which implements a bilinear interpolation based on assigning the weights of the pixels to their upper left coordinate corners
## alternatively, with the argument centroid = TRUE bilinear interpolation will be computed on the basis of centroids.
## depending on how data is passed in, Rcov requires different types of arguments to correctly harmonize the data with other covariate representations

## Unlike some other covariate types, Rcov can only handle covariates of the same class at any given time


Rcov <- function(..., interpolator = c('bilin', 'bin', 'spatstat')) {
    rasters <- list2(...)
    interpolator <- match.arg(interpolator, interpolator)
    object_types <- length(unique(lapply(rasters, class)))
    if (object_types != 1L) stop('Rcov requires all arguments be of the same type.\n Try using Rcov multiple times for each type of argument.')

    rasters <- lapply(rasters, convert_raster)

    new_rcrd(list(rasters = rasters),
             interpolator = interpolator,
             class = c('Rcov', 'spatial_covariate'))
}

evaluate.Rcov <- function(object, locations, ...){
    stopifnot(is_coord(locations))

    x <- coordx(locations)
    y <- coordy(locations)

    rasters <- field(object, 'rasters')
    interpolator <- switch(attr(object, 'interpolator'),
                           'bilin' = function(r) interp.im(r, x,y, bilinear = TRUE),
                           'bin' = r[x,y],
                           'spatstat' = function(r) interp.im(r, x,y, bilinear = FALSE))

    do.call(c, lapply(rasters, interpolator))
}

format.Rcov <- function(r, ...){
    rep('Raster Object', vec_size(r))
}

vec_ptype_abbr.Rcov <- function(r, ...){
    'Rcov'
}


convert_raster <- function(object, ...){
    UseMethod('convert_raster')
}

convert_raster.im <- function(object, ...){
    return(object)
}

convert_raster.default <- function(object, ...){
    cl <- class(object)[[1]]
    message <- paste0('Objects of type ', cl, ' not currently supported.')
    stop(message)
}

