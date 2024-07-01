#' covariate placeholders should carry the relevant information regarding
#' the appropriate single level entity  information needed in model fitting (such as mgcv::gam)
#'  cases include an age, sex, geo-coordinate (x,y), or possibly higher order coordinates, (x,y,t, w) for extra w
#' @export
placeholder_value <- function(value) {
  structure(value, class = c('placeholder_value', 'numeric'))
}
print.placeholder_value <- function(value) {
  cat(str_glue('(({value}))\n\n'))
}




## the goal of these is to get mgcv to work in a friendly way
## ' min.covariate_placeholder <- function(pl, ...) {min(extract_coords(pl)$x, ...)}
## ' max.covariate_placeholder <- function(pl, ...) {max(extract_coords(pl)$x, ...)}

#'  @export
covariate_placeholder <- function(data, coords) {
  structure(rep(placeholder_value(1), nrow(coords)),
            class = c('covariate_placeholder', 'numeric'),
            coords = coords,
            data = data)
}

print.covariate_placeholder <- function(object){
  cat('A Covariate Placeholder\n')
}

`[.covariate_placeholder` <- function(object, ...){
  covariate_placeholder(extract_data(object), extract_coords(object)[...,])
}

#' @export
extract_data <- function(object) {
  attr(object, 'data')
}

#' @export
extract_coords <- function(object) {
  attr(object, 'coords')
}

#' @export
remap <- function(placeholder, new_coords) {
  covariate_placeholder(extract_data(placeholder), new_coords)
}


#' Raster Type Covariates

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

  new_rcrd(list(rasters = rasters),
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
convert_raster <- function(object, ...){
  UseMethod('convert_raster')
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





#' @param object
#' @export
evaluate <- function(object, ...){
 UseMethod('evaluate')
}

evaluate.numeric <- function(object, locations){
    n <- nrow(locations)
    rep(object, each = n)
}

plot.spatial_covariate <- function(cov, W, ...){
    plot(spatstat.geom::as.im(function(x,y) evaluate(cov, coord(x,y), ...), W))
}
as.im.spatial_covariate <- function(cov, W, ...){
        spatstat.geom::as.im(as_Fcov(function(x,y) evaluate(cov, coord(x,y), ...), W))
}

evaluate.spatial_covariate <- function(object, ...){
    stop('Evaluate method not implemented')
}
## @param ... <[`dynamic-dots`][rlang::dyn-dots]> What these dots do.


#
#' Preparation of Convolutional Covariate Representations
#'
#' @param object
#' @param W
#' @param dimyx
#' @param fractional
#' @param normalize
#'
#' @return
#' @export
#'
#' @examples
conv_prepare <- function(object, W, dimyx, fractional, normalize){
    object$window <- spatstat.geom::union.owin(object$window, W)

    ## DivideByPixelArea guarantees that the integral of the rasterized process equals the
    ## actual value of the process, e.g. perfom normalization
    immat <- spatstat.geom::pixellate(object, DivideByPixelArea = normalize, fractional = fractional)
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


    list(covariate = for_conv,
         half_dists = half_dists,
         window = object$window,
         dims = dimyx,
         distances = dists,
         angle = angle
         )
}



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
Pcov <- function(..., W , dimyx, fractional = FALSE ){
  covars <- rlang::list2(...)
  prepped <- vector('list', length = length(covars))
  for (i in seq_along(covars)){
    covar <- covars[[i]]

    prepped[[i]] <- Pcov_prepare(covar, W, dimyx, fractional)

  }

  new_rcrd(list(pcov = prepped),
           class = c('Pcov', 'spatial_covariate'))
}

format.Pcov <- function(object, ...){
  rep('Pcov', vec_size(object))
}

Pcov_prepare <- function(object, W, dimyx, fractional){
  UseMethod('Pcov_prepare')
}

## Pcov should be evaluable when provided with a parametrically chosen kernel
evaluate.Pcov <- function(object, locations, kernel, ...){

  evaluate_single <- function(p, kernel, locations){
    conved <- fft(fft(kernel(p$distances)) * p$covariate, inverse = TRUE) / prod(p$dims)
    pim <- as.im(Re(matrix(conved[1:p$dims[[1]], 1:p$dims[[2]]])), W = p$window, dimyx = p$dimyx)
    interp.im(pim, coordx(locations), coordy(locations))
  }

  pred <- map(field(object, 'pcov'), function(p) evaluate_single(p, kernel, locations))
  reduce(pred, c)
}


setup <- function(object, locations, ...){
  UseMethod('setup')
}

## setup should provide a method to construct the pseudodata necessary for non-parametric basis setup
## currently only working for the single covariate case
## this needs to be reworked by rewriting placeholders to use vctrs
setup.Pcov <- function(object, locations){
  covariate_placeholder(
    object,
    cbind(coordx(locations), coordy(locations)))
}

Pcov_prepare.ppp <- function(object, W, dimyx, fractional){
  conv_prepare(object, W, dimyx, fractional, normalize = TRUE)
}


Pcov_prepare.default <- function(object, ...){
  cl <- class(object)[[1]]
  message <- paste0('Objects of type ', cl, ' not currently supported.')
  stop(message)
}




##Fcov



## convert to Fcov
as_Fcov <- function(object, ...){
  UseMethod('as_Fcov')
}

as_Fcov.Fcov <- function(object, ...){
  object
}

as.im.Fcov <- function(object, W) {
  as.im(function(x,y) evaluate(object, coord(x,y)), W = W)
}
as_Fcov.Pcov <- function(object, kernel, ...){
  force(kernel)
  f <- function(x,y) evaluate(object, locations = coord(x,y), kernel = kernel)
  Fcov(f)
}


as_Fcov.Lcov <- function(object, kernel, ...){
  force(kernel)
  f <- function(x,y) evaluate(object, locations = coord(x,y), kernel = kernel)
  Fcov(f)
}




## Fcov is similar in spirit to Rcov, but handles data via a functional representation
## this may be under some circumstances more amenable to handling uncertainty due to measurement error
## Fcov requires a basis or model for the representation of the data
Fcov <- function(..., formula, method, family){
  object <- list2(...)
  fitted <- lapply(object, function(ob) Fcov_prepare(ob, ...))

  new_rcrd(list(fitted = fitted), class = c('Fcov', 'spatial_covariate'))
}

format.Fcov <- function(ob, ...){
  rep('Fcov', vec_size(ob))
}


`+.Fcov` <- function(ob, ob2) {
  if (inherits(ob2, 'numeric')) {
    result <- map2(field(ob, 'fitted'),
                   ob2,
                   function(f1, f2) function(x,y) f1(x,y) + f2)
    return(Fcov(!!!result))
  }
  if(inherits(ob2, 'Fcov')){

    result <- map2(field(ob, 'fitted'),
                   field(ob2, 'fitted'),
                   function(f1, f2) function(x,y) f1(x,y) + f2(x,y))
    return(Fcov(!!!result))
  }
  stop('Not implemented')
}


`*.Fcov` <- function(ob, ob2) {
  if (inherits(ob2, 'numeric')) {
    result <- map2(field(ob, 'fitted'),
                   ob2,
                   function(f1, f2) function(x,y) f1(x,y) * f2)
    return(Fcov(!!!result))
  }
  if(inherits(ob2, 'Fcov')){

    result <- map2(field(ob, 'fitted'),
                   field(ob2, 'fitted'),
                   function(f1, f2) function(x,y) f1(x,y) * f2(x,y))
    return(Fcov(!!!result))
  }
  stop('Not implemented')
}

`-.Fcov` <- function(ob, ob2) {
  if (inherits(ob2, 'numeric')) {
    result <- map2(field(ob, 'fitted'),
                   ob2,
                   function(f1, f2) function(x,y) f1(x,y) - f2)
    return(Fcov(!!!result))
  }
  if(inherits(ob2, 'Fcov')){

    result <- map2(field(ob, 'fitted'),
                   field(ob2, 'fitted'),
                   function(f1, f2) function(x,y) f1(x,y) - f2(x,y))
    return(Fcov(!!!result))
  }
  stop('Not implemented')
}

`/.Fcov` <- function(ob, ob2) {
  if (inherits(ob2, 'numeric')) {
    result <- map2(field(ob, 'fitted'),
                   ob2,
                   function(f1, f2) function(x,y) f1(x,y) / f2)
    return(Fcov(!!!result))
  }
  if(inherits(ob2, 'Fcov')){

    result <- map2(field(ob, 'fitted'),
                   field(ob2, 'fitted'),
                   function(f1, f2) function(x,y) f1(x,y) / f2(x,y))
    return(Fcov(!!!result))
  }
  stop('Not implemented')
}




evaluate.Fcov <- function(ob, locations){
  x <- coordx(locations)
  y <- coordy(locations)
  preds <- map(field(ob, 'fitted'), function(f) f(x,y))
  reduce(preds, c)
}

setup.Fcov <- evaluate.Fcov


Fcov_prepare <- function(object, ...){
  UseMethod('Fcov_prepare')
}

Fcov_prepare.function <- function(object, ...) {
  stopifnot(names(formals(object)) == c('x', 'y'))
  object
}
## basis_select should, given a particular basis type, be able to pick the relevant basis elements and perform any sort of low-rank computations required
## basis select needs to be able to consider selecting between various types of basis, (cubic spline, thin plate, wavelet, fpca, gaussian process)
## it also needs to be able to make basis constructions either jointly across all covariate values at the same time (pooling the data)
## it should also be able to make pooling decisions using random subsets of the pooled data (for representativeness while also reducing computation costs).
## Alternatively it should be able to select bases for each model individually
## A final alternative is to use a specified portion of the data to make a single global basis selection

Fcov_prepare.data.frame <- function(data, formula, method, family){
  method(formula, family = family, data = data)
}
basis_select <- function(...){

}


## basis_project should perform model fitting to get basis representations of each of the functional covariates
## It should be able to fit each portion entirely seperately with its own smoothing parameter selection
## Alternatively, it should be additionally possible to pool smoothing parameter selection
## Finally, it should also be possible to fit models with a single totally pooled model across all covariate data (such as a large spatiotemporal regression)



## evaluate should be able to evaluate the fitted functional components at a set of specified coordinates

## Arguments to Fcov could be reasonably provide as either a list of data frames corresponding to the component regressions,
## or a single data frame corresponding to the total regression

## there could potentially be different basis selections for different components, e.g. if data is highly non-stationary over time

## new_rcrd(fields = list(coefs = coefs,
##                        ## the actual basis representation, coupled with the ability to interpret the basis at new points is critical
##                        which_basis,
##                        ## could possibly need influence function for down the line work
##                        ## this may be better computed via autodiff or explicit derivatives
##                        influence = influence,
##                        ## keeping the data as is would be way too memory intensive
##                        ## may need to refer to data in order to be able to calculate confidence intervals down the line
##                        data = data))
## ## the attributes of the record should contain the basis, or bases if different bases are used at varying times


##Ecov
## Ecov is for representing exposure based covariates.
## These are often used in general in ways that are quite similar to Fcovs, but may be more appropriate for modeling as the outcome


#' @export
Ecov <- function(..., outcome) {
  object <- rlang::list2(...)
  fitted <- lapply(object, function(ob) Ecov_prepare(ob, ...))

  vctrs::new_rcrd(list(outcome), class = c('Ecov', 'spatial_covariate'))
}

#' @export
format.Ecov <- function(ob, ...){
  rep('Ecov', vctrs::vec_size(ob))
}




##Lcov


## Lcov is meant for processing of linear network data
## an object of type 'linnet' or sf MULTILINE objects can be handled as data input
## alternatively, a collection of LINE or MULTILINE objects can also be handled as data input

## Lcov can be used to process multiple covariates at the same time, each covariate will be given the name of the argument if provided,
## or otherwise will simply be assigned the name of the symbol passed in
## Lcov returns a representation of the point process containing the necessary components for setting up a convolutional basis representation.
## The return value is a dummy vector of length 1 containing only placeholder numeric data (this may be used later to store genuinely useful information).
## The attributes contain angle orientation and distance matrices for setup of the convolutional basis. It also contains the precomputed fft of the
## pixellated version of the process
## additional attributes include the oservation window of the process and the dimensions of the process

Lcov <- function(..., W = NULL, dimyx = c(128, 128)){
  covars <- rlang::list2(...)

  for (i in seq_along(covars)){
    covar <- covars[i]

    covar[[i]] <- Lcov_prepare(covar[[i]], W, dimyx, fractional)
  }

  new_rcrd(list(pcov = covars),
           class = c('Lcov', 'spatial_covariate'))
}

## Lcov should be evaluable when provided with a parametrically chosen kernel
evaluate.Lcov <- function(object, locations, kernel, ...){

}


format.Lcov <- function(object, ...){
  rep('Lcov', vec_size(object))
}


Lcov_prepare.default <- function(object, ...){
  cl <- class(object)[[1]]
  message <- paste0('Objects of type ', cl, ' not currently supported.')
  stop(message)
}


Lcov_prepare <- function(object, W, dimyx, fractional){
  UseMethod('Lcov_prepare')
}



Lcov_prepare.linnet <- function(object, W, dimyx, fractional){
  conv_prepare(object, W, dimyx, fractional, normalize = FALSE)
}
