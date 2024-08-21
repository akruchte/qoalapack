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
Pcov <- function(..., dimyx = c(128, 128), W = NULL) {
    prepped <- lapply(list2(...), \(ob) pcov_prepare(ob, dimyx = dimyx, W = W))
    vctrs::new_vctr(prepped,
                    class = c('Pcov', 'spatial_covariate'))
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
    lapply(object, \(ob) prepare(ob, ...))
}


#' @exportS3Method
pcov_prepare.ppp <- function(object, dimyx, fractional, W, ...){
      conv_prepare(object, dimyx = dimyx, W = W)
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
conv_prepare <- function(object, dimyx, W){

    if (is.null(W)) W <-get_window_from_object(object)
    
    fractional = TRUE
    normalize = TRUE
    ## W as a possible buffer region?
    ## object$window <- spatstat.geom::union.owin(object$window, W)

    ## DivideByPixelArea guarantees that the integral of the rasterized process equals the
    ## actual value of the process, e.g. perfom normalization
    immat <- spatstat.geom::pixellate(object, DivideByPixelArea = TRUE, fractional = TRUE, dimyx = dimyx)
    ## establish coordinates of convolution raster
    ## the convolution raster extends out twice past the window in order to prevent circular convolution
    xcoords <- seq(from = immat$xrange[[1]], to = 2 * immat$xrange[[2]], length.out = 2 * dimyx[[1]])
    ycoords <- seq(from = immat$yrange[[1]], to = 2 * immat$yrange[[2]], length.out = 2 * dimyx[[2]])


    ## TODO alternative non fft based implementation

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
         angle = angle,
         unique_dists = unique(dists)
         )
}

 
get_max_distance <- function(ob) {
    stopifnot(inherits(ob, 'Pcov'))

    max_dist <- -Inf
    for (el in ob) {
        max_dist <- max(max_dist, max(el$unique_dists))
    }
    max_dist
}
