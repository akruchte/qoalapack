## library(rlang)
## library(vctrs)
## library(spatstat.core)
## source('coord.R')
## source('covariate_types.R')



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
    covars <- list2(...)

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



