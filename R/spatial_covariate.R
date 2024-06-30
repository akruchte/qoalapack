## covariate representations should have methods that make them easily coercible for representation in modeling
## specifically, there should be a method for evaluating the covariate representation at the points of a quadrature scheme


## source('coord.R')

#' @export
evaluate <- function(object, ...){
 UseMethod('evaluate')
}

#' @export
plot.spatial_covariate <- function(cov, W, ...){
    plot(as.im(function(x,y) evaluate(cov, coord(x,y), ...), W))
}

#' @export
as.im.spatial_covariate <- function(cov, W, ...){
        as.im(as_Fcov(function(x,y) evaluate(cov, coord(x,y), ...), W))
}

#' @export
evaluate.spatial_covariate <- function(object, ...){
    stop('Evaluate method not implemented')
}
## @param ... <[`dynamic-dots`][rlang::dyn-dots]> What these dots do.


## conv_prepare is a utility function used for preparing convolutional representations.
## it is used for both Lcov and Pcov
## Lcov and Pcov are essentially identical in representation and share most of the same code, but they are sufficiently
## different in interpretation and use that they are divided into completely seperate classes
#' @export
conv_prepare <- function(object, W, dimyx, fractional, normalize, .direction.preferred = NULL){
    if (!is.null(.direction.preferred)) stop("Preferred direction not implemented")

    object$window <- union.owin(object$window, W)

    ## DivideByPixelArea guarantees that the integral of the rasterized process equals the
    ## actual value of the process, e.g. perfom normalization
    immat <- pixellate(object, DivideByPixelArea = normalize, fractional = fractional)
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

    dists <- fft2shift(dists)
    ## setup angle orientation matrix
    angle <- outer(xcoords - median(xcoords),
                   ycoords - median(ycoords),
                   atan2)
    angle <- fft2shift(angle)


    list(covariate = for_conv,
         window = object$window,
         dims = dimyx,
         distances = dists,
         angle = angle
         )
}




