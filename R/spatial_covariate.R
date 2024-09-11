#' Spatial Covariate Classes
#'
#' @name spatial_covar



#' @export
evaluate.spatial_covariate <- function(object, ...){
    stop('Evaluate method not implemented')
}
## @param ... <[`dynamic-dots`][rlang::dyn-dots]> What these dots do.


placeholder_recycle <- function(object, target) {

}


#' @export
plot.spatial_covariate <- function(cov, W, ...){
    plot(spatstat.geom::as.im(function(x,y) evaluate(cov, coord(x,y), ...), W))
}

#' @export
as.im.spatial_covariate <- function(cov, W, ...){
        spatstat.geom::as.im(as_Fcov(function(x,y) evaluate(cov, coord(x,y), ...), W))
}


#' @exportS3Method
format.spatial_covariate <- function(object, ...){
    type <- class(object)[[1]]
    rep(type, vctrs::vec_size(object))
}


strip_meta <- c


default_spatial_points_config <- list(npoints = 100L, F2 = 10)


#' Generate Default Evaluation Points
#'
#' This function generates a set of evaluation points at a default set of locations.
#' @param An object of type extent or with an extent method
default_points <-function(ext) {
    SobolSequence::sobolSequence.points(
                       2,
                       default_spatial_points_config$F2,
                       default_spatial_points_config$npoints
                   )
}








##########
## PCOV ##
##########

#' Pcov
#' Pcov is meant for point process valued data,
#' an object of type 'ppp' or sf POINT objects can be handled as data input

#' Pcov can be used to process multiple covariates at the same time, each covariate will be given the name of the argument if provided,
#' or otherwise will simply be assigned the name of the symbol passed in
#' Pcov returns a representation of the point process containing the necessary components for setting up a convolutional basis representation.
#' The return value is a dummy vector of length 1 containing only placeholder numeric data (this may be used later to store genuinely useful information).
#' The attributes contain angle orientation and distance matrices for setup of the convolutional basis. It also contains the precomputed fft of the
#' pixellated version of the process
#' additional attributes include the observation window of the process and the dimensions of the process
#' a number of vctrs methods need to still be added for Pcovs, such as restore methods etc
#' @export
## W and dimyx should be moved into the attributes of the vector
## likewise distance and angle information should be shared between all covariates
Pcov <- function(..., resolution = c(128, 128), extent = NULL) {
    prepped <- lapply(list2(...), \(ob) pcov_prepare(ob, resolution, extent))
    vctrs::new_vctr(prepped,
                    class = c('Pcov', 'spatial_covariate'))
}

Pcov_impl <- function (points, extent) {
    
}

## Pcov should be evaluable when provided with a parametrically chosen kernel
#' @exportS3Method
evaluate.Pcov <- function(object, locations, kernel, ...){

    if (length(x) == 1) {
        p <- object[[1]]
        conved <- fft(fft(kernel(p$distances)) * p$covariate, inverse = TRUE) / prod(p$dims)
        pim <- spatstat.geom::as.im(Re(matrix(conved[1:p$dims[[1]], 1:p$dims[[2]]])), W = p$window, dimyx = p$dimyx)
        out <- interp.im(pim, coordx(locations), coordy(locations))
    }
    else {
           pred <- lapply(field(object, 'pcov'), function(p) evaluate(p, kernel, locations))
           out <- do.call(c, pred)
    }
    return (out)
   
}

#' Prepare
#'
#' generic interface for post initialization preparation of spatial covariates
#' @export
pcov_prepare <- function(object,  ...){
  UseMethod('pcov_prepare')
}



#' @exportS3Method
pcov_prepare.list <- function(object, ...){
    lapply(object, \(ob) pcov_prepare(ob, ...))
}


#' @exportS3Method
pcov_prepare.ppp <- function(object, resolution, extent, ...){
      conv_prepare(object, resolution, extent)
}

#' @exportS3Method
pcov_prepare.default <- function(object, ...){
  cl <- class(object)[[1]]
  message <- paste0('Objects of type ', cl, ' not currently supported.')
  stop(message)
}





get_window_from_object <- function(object)
{
    object$window
}


    
#' Preparation of Convolutional Covariate Representations
#' @export
conv_prepare <- function(object, resolution, extent, fractional = TRUE, normalize = TRUE){

    if (is.null(extent)) extent <-get_window_from_object(object)

    dimyx <- resolution
    ## W as a possible buffer region?
    ## object$window <- spatstat.geom::union.owin(object$window, W)

    ## DivideByPixelArea guarantees that the integral of the rasterized process equals the
    ## actual value of the process, e.g. perfom normalization
    immat <- spatstat.geom::pixellate(object, DivideByPixelArea = TRUE, fractional = TRUE, dimyx = dimyx)

    list(covariate =matrix(immat$v, nrow = resolution[1], ncol = resolution[2]),
         window = object$window,
         dims = dimyx
         )
}




##########
## RCOV ##
##########

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

    x <- coordx(locations)$x
    y <- coordy(locations)$y

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




##########
## LCOV ##
##########


#' Linear Network Covariates

#' Lcov is meant for processing of linear network data
#' an object of type 'linnet' or sf MULTILINE objects can be handled as data input
#' alternatively, a collection of LINE or MULTILINE objects can also be handled as data input

#' Lcov can be used to process multiple covariates at the same time, each covariate will be given the name of the argument if provided,
#' or otherwise will simply be assigned the name of the symbol passed in
#' Lcov returns a representation of the point process containing the necessary components for setting up a convolutional basis representation.
#' The return value is a dummy vector of length 1 containing only placeholder numeric data (this may be used later to store genuinely useful information).
#' The attributes contain angle orientation and distance matrices for setup of the convolutional basis. It also contains the precomputed fft of the
#' pixellated version of the process
#' additional attributes include the oservation window of the process and the dimensions of the process
#' @export
Lcov <- function(..., dimyx = c(128, 128)){
  
  covars <- rlang::list2(...)

  for (i in seq_along(covars)){
    covar <- covars[i]

    covar[[i]] <- Lcov_prepare(covar[[i]], W, dimyx)
  }

  vctrs::new_rcrd(list(pcov = covars),
           class = c('Lcov', 'spatial_covariate'))
}

## Lcov should be evaluable when provided with a parametrically chosen kernel
#' @export
evaluate.Lcov <- function(object, locations, kernel, ...){

}

#' @export
format.Lcov <- function(object, ...){
  rep('Lcov', vec_size(object))
}

#' @export
Lcov_prepare.default <- function(object, ...){
  cl <- class(object)[[1]]
  message <- paste0('Objects of type ', cl, ' not currently supported.')
  stop(message)
}

#' @export
Lcov_prepare <- function(object, W, dimyx){
  UseMethod('Lcov_prepare')
}


#' @export
Lcov_prepare.linnet <- function(object, W, dimyx){
  conv_prepare(object, dimyx)
}

#' @export
as_Lcov.im <- function(immat) {
    dimyx <- immat$dim
    ## establish coordinates of convolution raster
    ## the convolution raster extends out twice past the window in order to prevent circular convolution
    xcoords <- seq(from = immat$xrange[[1]], to = 2 * immat$xrange[[2]], length.out = 2 * dimyx[[1]])
    ycoords <- seq(from = immat$yrange[[1]], to = 2 * immat$yrange[[2]], length.out = 2 * dimyx[[2]])


    ## construct matrix, zero pad, and then perform fft
    for_conv <- matrix(0, nrow = dimyx[[1]] * 2, ncol = dimyx[[2]] * 2)
    for_conv[1:dimyx[[1]], 1:dimyx[[2]]] <- immat$v
    for_conv <- fft(for_conv)


    ## fft2shift is used after each to prevent phase issues down the line
    ## performing the fft shift during construction guarantees correctness of convolutions later without requiring complex tracking of indices
    ## setup distances matrix
    dists <- outer(xcoords - median(xcoords),
                   ycoords - median(ycoords),
                   \(x,y) sqrt(x^2 + y^2))


    half_dists <- seq(from = 0, to = max(dists), length.out = ceiling(dim(dists)[[1]]/2))

    dists <- fft2shift(dists)
    ## setup angle orientation matrix
    angle <- outer(xcoords - median(xcoords),
                   ycoords - median(ycoords),
                   atan2)
    angle <- fft2shift(angle)

    window <- owin(immat$xrange, immat$yrange)

   out <- list(covariate = for_conv,
         half_dists = half_dists,
         window = window,
         dims = dimyx,
         distances = dists,
         angle = angle,
         unique_dists = unique(dists)
         )


    out <- vctrs::new_rcrd(list(pcov = list(out)),
                    class = c('Lcov', 'spatial_covariate'))
    out
}
