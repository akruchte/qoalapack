#' Model Context
#'
#' Model context is an R6 class for capturing intermediate results of running mgcv functions, and not repeating them unnecessarily.
#' @export
model_context <- R6Class('model_context',
                         public = list(
                             context = NA,
                             set_context = function(new_context) self$context <- new_context))


model_configuration <- R6Class('model_configuration', list())

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



#'
#' for_conv has been previously fourier transformed
#' This function carries out the logic of convolving, and restructuring the data
#' @export
convolve_basis <- function(basis, pp_covariate) {

    resolution <- pp_covariate$dims

    covar_ix <- 1:resolution[1]
    covar_iy <- 1:resolution[2]

    basis_ix <- 1:resolution[1]
    basis_iy <- 1:resolution[2]
    
    out_ix <- (resolution[1] %/% 2) + 1:resolution[1]
    out_iy <- (resolution[2] %/% 2) + 1:resolution[2]
    
    outbuf <- matrix(0, nrow = resolution[1], ncol = resolution[2])
    buf1 <- matrix(0, nrow = resolution[1] * 2L, ncol = resolution[2] * 2L)
    buf2 <- matrix(0, nrow = resolution[1] * 2L, ncol = resolution[2] * 2L)

    dims <- pp_covariate$dims
    window <- pp_covariate$window

    buf1[covar_ix, covar_iy] <- pp_covariate$covariate
    buf2[basis_ix, basis_iy] <- basis

    vec <- Re(fft(fft(buf1) * fft(buf2), inverse = TRUE)) / prod(2 * resolution)
    outbuf[] <- vec[out_ix, out_iy]
    Rcov(spatstat.geom::as.im(outbuf, W = window))
}



## required mgcv function
## this function expects to receive data in the form of a covariate placeholder.
## The levels are included directly in the data
## actual coordinates and pcov data is encoded internally in attributes


#' extract_basis_from_spec
#'
#' mgcv smooth specs take provided bases and store this information in the object
#' class appended with ".smooth.spec"
#' this simply reverses this
extract_basis_from_spec <- function(smooth_spec) {
    stringr::str_remove(class(smooth_spec), stringr::fixed(".smooth.spec"))
}


#' @exportS3Method
smooth.construct.conv.smooth.spec <- function(object, data, knots) {

    ldata <- data[[object$term]]
    
    extra <- object$xt
    ctxt <- extra$context
    

    term <- object$term
    
    coords <- extract_coords(ldata)
    conv_data <- extract_data(ldata)

    ## get points from process
    ## compute distance matrix from all points of process to all target points

    ## construct spline basis on this

    
    ## use the below logic for linear process covariates

    resolution <- conv_data[[1]]$dims

    xseq <- seq(from = -1L, to = 1L, length.out = resolution[1])
    yseq <- seq(from = -1L, to = 1L, length.out = resolution[2])



    distrast <- outer(xseq, yseq, function(x,y) sqrt(x^2 + y^2))
    local_data <- list(distances = c(distrast))
    
    nknots <- object$bs.dim
    if(nknots < 1) {
        nknots <- 10
        object$bs.dim <- nknots
    }
    
    n_right_boundary_knots <- 2L
    n_left_boundary_knots <- 2L
    
    lknots <- seq(from = 0L, to = 1L, length.out = (nknots + n_right_boundary_knots))
    interval <- lknots[2] - lknots[1] 
    lknots <- c(-2 * interval, -1 * interval, lknots)

    ## TODO user supplied knots
    knots <- list(distances = lknots)

    basis_term <- extract_basis_from_spec(object)[-1]
    if (is_empty(basis_term)) basis_term <- 'bs2'

    intern_call <- s(distances,  bs = basis_term, fx = object$fixed, k = object$bs.dim, xt = object$xt)
    intern_call$label <- paste0('conv(', object$term, ')')



    
    basis <- smooth.construct(intern_call, data = local_data, knots = knots)

    ## preallocate output design
    distance_design <- basis$X
    ## number of distinct Pcov values
    bases <- vector(mode = 'list', length = length(conv_data))    

    for (i in seq_along(bases)){
        current_pp <- conv_data[[i]]

        ncols <- ncol(distance_design)
        interp_basis <- vector(mode = 'list', length = ncols)

        for(basis_index in 1:ncols)
        {
            interp_basis[[basis_index]] <- convolve_basis (distance_design[,basis_index], current_pp)
        }
        bases[[i]] <- interp_basis
    }

    ## then in this step apply this to each marked set seperately
    ## in that way everything is now pooled

    object$internal_basis <- basis
    object$interpolation_basis <- bases
    class(object) <- 'Convspline.smooth'


    object$X <- Predict.matrix.Convspline.smooth(object, data)

    return(object)
}



## adaptive convolutions use a three dimensional convolution
## and then add an additional adaptive surface penalty

## required mgcv function

#' @export
Predict.matrix.Convspline.smooth <- function(object, data) {
    coords <- extract_coords(data[[object$term]])
    
    interp_basis <- object$interpolation_basis

    ncoord <- length(coords)
    
    nr <- ncoord * length(interp_basis)
    ## assumes all bases have same dimension which they certainly should
    nc <- length(interp_basis[[1]])
    Xmat <- matrix(0, nrow = nr, ncol = nc)

    for (i in 1:length(interp_basis)){

        locinterp <- interp_basis[[i]]

        for (j in 1:nc){
            ## e <<- environment()
            row_range <- (((i - 1) * ncoord) + 1) : (((i) * ncoord) )
            interp <- evaluate(locinterp[[j]], coords)
            Xmat[row_range, j] <- interp
        }
    }

    Xmat
}







alg_environment <- function ( ){
    e <-  rlang::env()

    resolution <- c(128, 128)

    e$cbuf_covariate <- matrix(0, nrow = resolution[1] * 2L, ncol = resolution[2] * 2L)
    e$cbuf_basis <- matrix(0, nrow = resolution[1] * 2L, ncol = resolution[2] * 2L)
    e$outbuf <- matrix(0, nrow = resolution[1], ncol = resolution[2])
    e$drast <- outer(
        X = seq(from = -1, to = 1, length.out = resolution[1]),
        Y = seq(from = -1, to = 1, length.out = resolution[2]),
        \(x,y) sqrt(x^2 + y^2))

    covar_ix <- 1:resolution[1]
    covar_iy <- 1:resolution[2]

    basis_ix <- 1:resolution[1]
    basis_iy <- 1:resolution[2]
    
    out_ix <- (resolution[1] %/% 2) + 1:resolution[1]
    out_iy <- (resolution[2] %/% 2) + 1:resolution[2]
}


#' sgam is a thin wrapper around gam that provides revision capabilities. 
## 
#' @export
sgam <- function(formula, data, conv_control, ... ) {

    resolution <- c(128, 128)

    cbuf_covariate <- matrix(0, nrow = resolution[1] * 2L, ncol = resolution[2] * 2L)
    cbuf_basis <- matrix(0, nrow = resolution[1] * 2L, ncol = resolution[2] * 2L)
    outbuf <- matrix(0, nrow = resolution[1], ncol = resolution[2])
    covar.test <- spatstat.geom::as.im(swedishpines)$v


    basis <- splines::spline.des(
                          knots = seq(from = -0.1, to = 2, length.out = 45),
                          c(drast),
                          outer.ok = TRUE)
                          
                                 


    cbuf_covariate[covar_ix, covar_iy] <- c(covar.test)
    ## line for visibility
    ## cbuf_covariate[12, covar_iy] <- 1
    cbuf_basis[basis_ix, basis_iy] <- c(basis$design[,28])

    x <- fft(cbuf_covariate)
    y <- fft(cbuf_basis)

    out <- Re(fft(x * y, inverse = TRUE)) / prod(2 * resolution)


    outbuf[] <- out[out_ix, out_iy]
    
    fit <- gam(formula, data, conv_control, ...)
    class(fit) <- c('sgam', class(fit))
}





## given a certain parametric model use a spline based model to test correctness of the specification
## update the parametric model in one direction or another on the basis of this test

## Once we have a correctly specified test according to the spline approximation we can then revise the spline model to look at max distance specifications or whatnot

spline_spec_test <- function ( ) {
    
}


## TODO
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

