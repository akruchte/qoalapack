## source('covariate_types.R')
## source('utilities.R')
## examples for documentation


test_outcome_mod <- function() {
    library(tidyverse)
    library(spatstat)
    set.seed(0)
    gen_covar <- function(...) spatstat.random::rpoispp(35)
    Ys <- map(1:5, gen_covar)
    Xs <- map(1:5, gen_covar)

    points <- reduce(Ys, superimpose)

    pcovs <- map(Xs, Pcov)

    pc <- covariate_placeholder(pcovs, coords(points))

    Yout <-     rnorm(163)

    df <- tibble(out = rnorm(815), cov = pc)

    gam(out ~ s(cov, bs = 'conv'), data = df)
}


## library(mgcv)
## library(pracma)

## library(spatstat)
## library(dplyr)
## library(purrr)
## library(stringr)
## library(glue)
## library(rlang)


## crude temporary implementation of outcome model
## takes in a fitted model and provides a simple wrapper for calculating
## quantities needed for estimation
#' @export
outcome_model <- function(model, prediction_locations) {
    list(model = model,
         mufun = outcome_fun(model, prediction_locations))
}

## returns a function f(i) for treatment A indexed by i.
## f evaluates the predicted intensity at each of the originally specified locations when assigned a treatment value of A
## expects that the model have treatment provided as the first term
#' @export
outcome_fun <- function(omodel, prediction_locations){
    if(missing(prediction_locations)) stop('Prediction dataset required')
    linkinv <- exp

    mu <- predict(omodel, newdata = prediction_locations, type = 'terms')

    treat_term <- mu[,1]

    mu_term <- rowSums(as.matrix(mu[,-1])) + attr(mu, 'constant')

    denominator_mu <- linkinv(mu)
    ## homoskedastic normal case

    ## mufun is a function that predicts the potential outcome, at covariate values corresponding to those at each and every of the provided prediction locations
    ## For each prediction location it returns the predicted value corresponding to the ith level of the observed treatments if the argument i is provide.
    ## if a is provided, it predicts on the basis of the value of a
    #' @export
    mufun <- function(i, a) {
        if(missing(i) & missing(a)){
            stop('Conditional mean prediction requires either an index referring to an observed treatment (i), or a specific value of treatment (a)')
        }
        if(!missing(i)){
            stopifnot(missing(a))
            return(exp(mu_term + treat_term[i]))
        }

        if(!missing(a)){
            stopifnot(missing(i))
            nd <- mutate(prediction_locations, treatment = a)
            predict(omodel, newdata = nd, type = 'response')
        }
    }

    return(mufun)

}

## formula extractor
## deconstruct formula and determine relevant terms, etc

#' @export
ordinal_model <- function(){
    gam(outcome ~ s(x,y), family = ocat(R = n))
}

## this should be the principle function
## it should receive datastructures representing the data, model specifications, control specs, prediction locations, and everything else
## needed for ultimate use of the outcome model in causal estimation
## it should also handle basis setup logic, and any other logic needed to make future use of the model proceed smoothly

## we previously had a function here called outcome data
## this should instead be a tibble representation of outcomes, and covariates

## outcome_control is used for specifying model specific and parameter estimation details

## quadrature_control can be a quadrature object, a set of specified NEED TO DETERMINE
## optimizer_control is used for selecting the optimization parameters, and method of optimization
## if method = 'gam' mgcv is used directly for estimation using the mgcv native defaults. Alternative methods may be provided
## down the line if distributed optimization is required
#' @export
outcome_control <- function(quadrature_control,
                            optimizer_control,
                            gam_control,
                            raster_control
                            ){

}

#' @export
ppmod <- function(Y, Q, ppcov, covariates = NULL, dimyx = c(128, 128), k = NULL, bs = 'tp') {
    stopifnot(!is.null(names(ppcov)))
    stopifnot(is.list(ppcov))

    prep <- mpl_prepare(Y, Q, ppcov, covariates, dimyx)
    ppcov_list <- names(ppcov)


    ppcov <- map_chr(ppcov_list, function(str) paste0('s(', str, ', bs = c("conv", bs), k = ', k,')'))

    cov_list <- names(covariates)
    cov_list <- if(!is.null(cov_list)) paste0(cov_list, collapse = '+')

    gam_formula <- paste0('outcome ~ ',  paste0(ppcov, collapse = ' + '))


    gam_formula <- paste0(c(gam_formula, cov_list), collapse = '+')

    fit <- suppressWarnings(gam(formula(gam_formula), data = prep, family = poisson(), weights = w))

    fit$Y <- Y
    fit$Q <- Q
    fit$gam_data <- prep
    fit$ppcovs <- ppcov_list
    fit$dim <- dimyx
    class(fit) <- c('ppmod', class(fit))
    fit

}

#' @export
print.ppmod <- function(object) {
    cat('A Point Process Model:\n\n')
    cat(glue('Outcome process with {npoints(object$Y)} points.\n\n'))
    cat(glue('Quadrature scheme with {npoints(object$Q$dummy)} control points.\n\n'))

    cat(glue('Point process valued covariates: {object$ppcov}.\n\n'))

}

#' @export
print.counterfactual <- function(object){
    cat('Counterfactual Modified Point Process: \n\n')
    print.ppmod(object)

    cat('\n\n')
    cat(glue('Counterfactually modified covariates: {object$counterfactual_covs}.\n\n'))
}

#' @export
predict.ppmod <- function(object, newdata, ... ) {

    if(missing(newdata)){
        return(predict.gam(object, ...))
    }

    for (covariate in object$ppcovs) {
        newdata[,covariate] <- remap(object$gam_data[,covariate],
                                     new_coords = newdata[,c('x', 'y')])
    }


    predict.gam(object, newdata = newdata, ...)
}

#' @export
prep_outcome <- function(Y, Q){
    bind_cols(
        bind_rows(
            bind_cols(coords(Q$data), outcome = 1),
            bind_cols(coords(Q$dummy), outcome = 0),
            ),
        w = Q$w) |>
        mutate(outcome = outcome / w)
}



## given point process and quadscheme setup for point process models
resp_value_prepare <- function(Y, Q) {
    prepped_data <- bind_cols(
        bind_rows(
            bind_cols(coords(Q$data), outcome = 1),
            bind_cols(coords(Q$dummy), outcome = 0),
            ),
        w = Q$w) |>
        mutate(outcome = outcome / w)

    prepped_data
}


## prepare data
## mpl prepare is now obsolete
## should be replaced with data frame based setup functions
## ppcov is a list of lists of point processes

#' @export
mpl_prepare <- function(Y, Q,  ppcov = NULL, covariates = NULL, dimyx = c(128, 128))
{
    prepped_data <- resp_value_prepare(Y, Q)
    
    prepped_pp_exposures <- vector(mode = 'list', length(ppcov))
    names(prepped_pp_exposures) <- names(ppcov)

    for(i in seq_along(ppcov)) {

        ## W = Y$window,
        imm <- Pcov(ppcov[[i]],  dimyx = dimyx)[[1]]
        res <- covariate_placeholder(imm, coords(Q))

        prepped_pp_exposures[[i]] <- res
    }

    prepped_covariates <- vector(mode = 'list', length(ppcov))
    names(prepped_covariates) <- names(covariates)
    
    for (i in seq_along(covariates)){
        covar <- covariates[[i]]
        if (is.im(covar)){
            ## new_covar <- interp.im(covar, x = prepped_data[,1], y = prepped_data[,2])
            new_covar <- covar[list(x = prepped_data[,1], y = prepped_data[,2])]
        }
        else if (is.function(covar)) {
            if (! identical(names(formals(covar)), c('x', 'y'))) stop('Functional covariates must have formal arguments x and y')
            new_covar <- covar(prepped_data[,1], prepped_data[,2])
        }

        else {
            stop('Covariates must be passed as either functions of x and y arguments or images')
        }
        prepped_covariates[[i]] <- new_covar

    }

    gam_data <- bind_cols(
        prepped_data,
        prepped_pp_exposures,
        prepped_covariates
    )

    
    gam_data
}

#' @export
update_exposure <- function(model, new_exposure) {

    newdata <- model$gam_data

    ## data checks
    stopifnot(is.list(new_exposure))
    if(any(is.null(names(new_exposure)))) stop('New exposures must have names')
    stopifnot(all(names(new_exposure) %in% model$ppcovs))


    which_covs <- model$ppcovs %in% names(new_exposure)
    covs <- model$ppcovs[which_covs]


    for (cov in covs) {
        newdata[[cov]] <- covariate_placeholder(
            Pcov(new_exposure[[cov]], model$Y, model$dim),
            extract_coords(model$gam_data[[cov]])
        )
    }
    smooths <- fit$smooth

    for(i in seq_along(smooths)){
        if(is(smooths[[i]], 'Convspline.smooth')){
            object <- smooths[[i]]
            cov <- object$term
            conv_data <- extract_data(newdata[[cov]])

            covlen <- length(conv_data)

            

            basis <- list( c(field(conv_data, 'pcov')[[1]]$distances))
            names(basis) <- object$term
            
            new_basis <- Predict.matrix(object$internal_basis, basis)

            
            object$interpolation_basis <- apply(
                new_basis, 2,
                function(basis) convolve_basis(basis,
                                               field(conv_data, 'pcov')[[1]]$covariate,
                                               field(conv_data, 'pcov')[[1]]$dims,
                                               field(conv_data, 'pcov')[[1]]$window,
                                               coords)
            )

            smooths[[i]] <- object
        }
    }
    model$smooth <- smooths
    model$counterfactual_covs <- covs

    if(!is(model, 'counterfactual')) class(model) <- c('counterfactual', class(model))

    model
}


## required mgcv function
## this function expects to receive data in the form of a covariate placeholder.
## The levels are included directly in the data
## actual coordinates and pcov data is encoded internally in attributes

#' @export
smooth.construct.conv.smooth.spec <- function(object, data, knots) {
    
    conv_data <- purrr::reduce(extract_data(data[[object$term]]), c)
    coords <- extract_coords(data[[object$term]])

    max_dist_prop <- object$xt$max_dist_prop
    if (is.null(max_dist_prop)) max_dist_prop <- 0.25


    term <- object$term

    ## if no secondary basis is provided default to p-splines
    basis_term <- 'ps'
        if (length(class(object)) > 1) {
            basis_term <- str_extract(class(object)[[2]], '[a-zA-Z]+')
        }

            
    intern_call <- s(distances,  bs = basis_term, fx = object$fixed, k = object$bs.dim)
    intern_call$label <- paste0('conv(', term, ')')

        ## local_data <- list(distances = unique(c(field(conv_data, 'pcov')[[1]]$distances)))
   ## TODO     make user configurable and provide better defaults
    distances <- seq(from = 0, to = max_dist_prop * max(unique(c(field(reduce(conv_data, c), 'pcov')[[1]]$distances))), length.out = 1000)
    local_data <- list(distances = distances)


        ## internal basis construction including penalty
    basis <- smooth.construct(intern_call, data = local_data, knots = knots)
    basis$og_data <- local_data


    ## I don't think this is necessary?
    pred_dat <- list(distances = reduce(map(field(reduce(conv_data, c), 'pcov'), 'distances'), c))


    basis$internal_basis <- basis
    basis$term <- object$term

    ## then in this step apply this to each marked set seperately
    ## in that way everything is now pooled

    local_conv <- function(to_conv, pcovar, coords) {
        covar <- pcovar$covariate
        dims <- pcovar$dims
        window <- pcovar$window
        convolve_basis(to_conv, covar, dims, window, coords)
    }
    
    bases <- vector(mode = 'list', length = length(conv_data))
    
    for (i in seq_along(bases)){
        lcd <- field(conv_data, 'pcov')[[i]]
        pred_dat <- list(distances = lcd$distances)
        Xloc <- Predict.matrix(basis$internal_basis, data = pred_dat)
        Xout <- 0 * Xloc

        ncols <- ncol(Xloc)
        interp_basis <- vector(mode = 'list', length = ncols)
        
        for(col in 1:ncols) {
            lc <- local_conv(Xloc[,col], lcd)

            interp_basis[[col]] <- lc
        }
        bases[[i]] <- interp_basis
        
    }

    ## then in this step apply this to each marked set seperately
    ## in that way everything is now pooled
    basis$interpolation_basis <- bases

    class(basis) <- 'Convspline.smooth'

    basis$X <- Predict.matrix.Convspline.smooth(basis, data)
    basis

}


smooth.construct.lconv.smooth.spec <- function(object, data, knots) {
    
    conv_data <- extract_data(data[[object$term]])
    coords <- extract_coords(data[[object$term]])

    max_dist_prop <- object$xt$max_dist_prop
    if (is.null(max_dist_prop)) max_dist_prop <- 0.25


    term <- object$term

    ## if no secondary basis is provided default to p-splines
    basis_term <- 'ps'
        if (length(class(object)) > 1) {
            basis_term <- str_extract(class(object)[[2]], '[a-zA-Z]+')
        }

            
    intern_call <- s(distances,  bs = basis_term, fx = object$fixed, k = object$bs.dim)
    intern_call$label <- paste0('conv(', term, ')')

        ## local_data <- list(distances = unique(c(field(conv_data, 'pcov')[[1]]$distances)))
   ## TODO     make user configurable and provide better defaults
    distances <- seq(from = 0, to = max_dist_prop * max(unique(c(field(reduce(conv_data, c), 'pcov')[[1]]$distances))), length.out = 1000)

    
    local_data <- list(distances = distances)


        ## internal basis construction including penalty
    basis <- smooth.construct(intern_call, data = local_data, knots = knots)
    basis$og_data <- local_data


    ## I don't think this is necessary?
    pred_dat <- list(distances = reduce(map(field(reduce(conv_data, c), 'pcov'), 'distances'), c))


    basis$internal_basis <- basis
    basis$term <- object$term

    ## then in this step apply this to each marked set seperately
    ## in that way everything is now pooled

    
        lcd <- field(conv_data, 'pcov')[[1]]
        pred_dat <- list(distances = lcd$distances)
        Xloc <- Predict.matrix(basis$internal_basis, data = pred_dat)
        Xout <- 0 * Xloc

        ncols <- ncol(Xloc)
        interp_basis <- vector(mode = 'list', length = ncols)

    
    for(j in 1:ncols) {

        covar <- lcd$covariate
        dims <- lcd$dims
        window <- lcd$window
        convolve_basis(Xloc[,j], covar, dims, window, coords)

            interp_basis[[col]] <- lc
        }
        
        


    ## then in this step apply this to each marked set seperately
    ## in that way everything is now pooled
    basis$interpolation_basis <- bases

    class(basis) <- 'Convspline.smooth'

    basis$X <- Predict.matrix.Convspline.smooth(basis, data)
    basis

}





## adaptive convolutions use a three dimensional convolution
## and then add an additional adaptive surface penalty



if(experiment <- FALSE){
    adapt_resolution <- 128
    max <-1
    adapt_scaling <-seq(from = 0.01, to = max, length.out = adapt_resolution)

    surf <- as.matrix(as.im(swedishpines))
    arrk <- array(data = 0, dim = c(128, 128, 128))
    adaptconv <- arrk
    kerns  <- arrk
    for (k in 1:adapt_resolution) {
        arrk[,,k] <-surf
        scale <-adapt_scaling[k]
        dist <- outer(seq(from = -1, to = 1, length.out = 128),
                      seq(from = -1, to = 1, length.out = 128),
                      function(x,y) exp(sqrt(scale * (x^2 + y^2))))

        kerns[,,k] <- dist
        ## problem here but it'll work for now
        adaptconv[,,k] <- Re(fft2shift(fft(fft(surf) * fft(dist), inverse = TRUE)))
        trilinear_interp <- function(x,y, z){
            ## get coords and do trilinear interpolation
        }

        ## z is a height function of time that determines convolutional scale.
        ## We extend the model with regularity terms capturing adaptive structure.


    }


    Y <- swedishpines
    Q <- quadscheme(Y)
    mpl_prepare(Y, Q, ppcov = Y, covariates = NULL)




    ## then do slice by slice convolution
    ## is there a reason to consider 3d-convolution? It might make sense.

    ## Then the resulting convolutional kernel can be identified by evaluating according to some varietal structure.
    ## We'll assume that K(D, \alpha) where alpha is itselfa function of space.


}


## required mgcv function

#' @export
Predict.matrix.Convspline.smooth <- function(object, data) {

    coords <- extract_coords(data[[object$term]])

    interp_basis <- object$interpolation_basis
    ncoord <- nrow(coords)
    nr <- ncoord * length(interp_basis)
    nc <- length(interp_basis[[1]])
    Xmat <- matrix(0, nrow = nr, ncol = nc)


    for (i in 1:length(interp_basis)){
        locinterp <- interp_basis[[i]]
        
        for (j in 1:nc){
            
            row_range <- (((i - 1) * ncoord) + 1) : (((i) * ncoord) )
            
            Xmat[row_range, j] <- spatstat.geom::interp.im(locinterp[[j]], coords)
        }
    }

    Xmat


    
}

#' @export
convolve_basis <- function(basis, for_conv, dims, window, coords){
    basis <- fft(basis)
    dim(basis) <- dims * 2


    vec <- Re(fft(basis * for_conv, inverse = TRUE) / prod(dim(basis)))
    spatstat.geom::as.im(vec[1:dims[1], 1:dims[2]], W = window)
}


