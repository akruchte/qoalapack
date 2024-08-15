## Raster Type Covariates

#'  Rcov handles data passed in raster formats.
#' Currently, only 'im' type objects are handled
#' Rcov does minimal handling of the raster objects, its primary goal is for metadata management and allowing natural manipulation of raster valued covariates

#' @param interpolator
#' options include 'bin' which treats the raster as a piecewise constant surface given by the raster values,
#' and 'bilin' which implements a bilinear interpolation based on assigning the weights of the pixels to their upper left coordinate corners
#' alternatively, with the argument centroid = TRUE bilinear interpolation will be computed on the basis of centroids.
##depending on how data is passed in, Rcov requires different types of arguments to correctly harmonize the data with other covariate representations




#' @export
Rcov <- function(..., interpolator = c('bilin', 'bin', 'spatstat')) {
  rasters <- rlang::list2(...)
  interpolator <- match.arg(interpolator, interpolator)
  object_types <- length(unique(lapply(rasters, class)))
  
  if (object_types != 1L) stop('Rcov requires all arguments be of the same type.\n Try using Rcov multiple times for each type of argument.')

  rasters <- lapply(rasters, convert_raster)

  vctrs::new_rcrd(list(rasters = rasters),
           interpolator = interpolator,
           class = c('Rcov', 'spatial_covariate'))
}

#' @export
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

#' @export
format.Rcov <- function(r, ...){
  rep('Raster Object', vec_size(r))
}

#' @export
vec_ptype_abbr.Rcov <- function(r, ...){
  'Rcov'
}


#' @export
convert_raster.im <- function(object, ...){
  return(object)
}

#' @export
convert_raster.default <- function(object, ...){
  cl <- class(object)[[1]]
  message <- paste0('Objects of type ', cl, ' not currently supported.')
  stop(message)
}

