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
Pcov <- function(..., dimyx = c(128, 128)) {
    prepped <- lapply(list2(...), \(ob) prepare(ob, dimyx = dimyx))
    vctrs::new_vctr(prepped,
                    class = c('Pcov', 'spatial_covariate'))
}



#' @exportS3Method
prepare.list <- function(object, ...){
    lapply(object, \(ob) prepare(ob, ...))
}


#' @exportS3Method
prepare.ppp <- function(object, dimyx, fractional){
      conv_prepare(object, dimyx = dimyx)
}

#' @exportS3Method
prepare.default <- function(object, ...){
  cl <- class(object)[[1]]
  message <- paste0('Objects of type ', cl, ' not currently supported.')
  stop(message)
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





