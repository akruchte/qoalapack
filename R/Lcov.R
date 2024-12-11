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
