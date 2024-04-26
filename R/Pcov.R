library(rlang)
library(vctrs)
library(spatstat.core)
source('utilities.R')
source('covariate_placeholder.R')
source('coord.R')
source('covariate_types.R')



## a number of vctrs methods need to still be added for Pcovs, such as restore methods etc

## Pcov is meant for point process valued data,
## an object of type 'ppp' or sf POINT objects can be handled as data input

## Pcov can be used to process multiple covariates at the same time, each covariate will be given the name of the argument if provided,
## or otherwise will simply be assigned the name of the symbol passed in
## Pcov returns a representation of the point process containing the necessary components for setting up a convolutional basis representation.
## The return value is a dummy vector of length 1 containing only placeholder numeric data (this may be used later to store genuinely useful information).
## The attributes contain angle orientation and distance matrices for setup of the convolutional basis. It also contains the precomputed fft of the
## pixellated version of the process
## additional attributes include the observation window of the process and the dimensions of the process


## W and dimyx should be moved into the attributes of the vector
## likewise distance and angle information should be shared between all covariates
Pcov <- function(..., W , dimyx, fractional = FALSE ){
    covars <- list2(...)
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

