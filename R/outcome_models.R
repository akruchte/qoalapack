#' @exportS3Method
print.Convspline.smooth <- function(ob, ...) {
    cat("A Convolutional Spline Basis: \n ")
    cat("Object names:\n")
    cat(paste(names(ob), collapse = "   "))
    cat("\n")
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
convolve_basis <- function(basis, pp_covariate, resolution) {

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
    odata <- data
    data <- data[[object$term]]
    
    extra <- object$xt
    max_distance_prop <- extra$max_distance_prop
    use_regularity <- extra$use_regularity
    regularity_model <- extra$regularity_model

    if(is.null(use_regularity)) use_regularity <- FALSE
     if(is.null(regularity_model)) {
        regularity_model <- "distance_gaussian"
    } 
    
    if(is.null(max_distance_prop)) max_distance_prop <- 1
    ctxt <- extra$context

    resolution <- object$xt$resolution
    if (is.null(resolution)) resolution <- c(128, 128)

    term <- object$term
    
    coords <- extract_coords(data)
    conv_data <- extract_data(data)


    if (regularity_model == "2d_gaussian") {
        xseq <- seq(from = -1L, to = 1L, length.out = resolution[1])
        yseq <- seq(from = -1L, to = 1L, length.out = resolution[2])

        xygrid <- expand_grid(y = yseq, x = xseq)
        distrast <- outer(xseq, yseq, function(x,y) sqrt(x^2 + y^2))
        clip <- distrast <= max_distance_prop
        local_data <- list(x = xygrid$x, y = xygrid$y)
        
        nknots <- object$bs.dim

        if(nknots < 1) {
            nknots <- 10
            object$bs.dim <- nknots
        }
        
        lknots <- seq(from = -max_distance_prop, to = max_distance_prop, length.out = nknots)

        ## TODO user supplied knots
        knots <- list(x = lknots, y = lknots)


        basis_term <- extract_basis_from_spec(object)[-1]
        if (is_empty(basis_term)) basis_term <- 'bs2'

        intern_call <- te(x,y,  fx = object$fixed, k = c(object$bs.dim, object$bs.dim))
        intern_call$label <- paste0('conv(', object$term, ')')

        basis <- smooth.construct(intern_call, data = local_data, knots = knots)

        ## preallocate output design
        distance_design <- basis$X
    }

    if (regularity_model == "distance_gaussian") {
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
        
        lknots <- seq(from = 0L, to = max_distance_prop * sqrt(2), length.out = (nknots + n_right_boundary_knots))
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

    }


    current_pp <- conv_data[[1]]

    ncols <- ncol(distance_design)
    interp_basis <- vector(mode = 'list', length = ncols)

    for(basis_index in 1:ncols)
    {
        interp_basis[[basis_index]] <- convolve_basis (distance_design[,basis_index], current_pp, resolution)
    }



    ## then in this step apply this to each marked set seperately
    ## in that way everything is now pooled

    object$internal_basis <- basis
    object$interpolation_basis <- interp_basis
    object$evaluation_coords <- coords
    object$window <- conv_data[[1]]$window
    object$im_dims <- conv_data[[1]]$dims
    object$local_data <- local_data
    object$use_regularity <- use_regularity
    object$regularity_model <- regularity_model
    
    class(object) <- 'Convspline.smooth'

    object$X <- Predict.matrix.Convspline.smooth(object, odata)

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
    
    nr <- ncoord 
    ## assumes all bases have same dimension which they certainly should
    nc <- length(interp_basis)
    Xmat <- matrix(0, nrow = nr, ncol = nc)

    for (i in 1:length(interp_basis)){

        locinterp <- interp_basis[[i]]

        interp <- spatstat.geom::interp.im(locinterp[[1]], cbind(coordx(coords), coordy(coords)))
        Xmat[, i] <- interp
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
updateF <- function(F, G,  offset){
    G$offset <- offset
    F <- gam(G = G, start = coef(F))
    class(F) <- c("indirect_gam", class(F))
    F
}


#' export
indirect_gam <- function (
                          formula,
                          data,
                          family,
                          indir.max.iter = 30,
                          delta.tol = 1e-6,
                          ...)
{
    gam_ <- gam
    ## G <- gam(formula = formula, family = family, data = data, offset = rep(0, nrow(data)), fit = FALSE)
    G <- gam_(formula = formula, family = family, data = data, offset = rep(0, nrow(data)), fit = FALSE, ...)
    F <- gam(G = G)
    class(F) <- c("indirect_gam", class(F))

    smooths <- F$smooth
    conv_smooths <- keep(smooths, \(s)inherits(s, "Convspline.smooth") )
    conv_smooths <- keep(conv_smooths, \(s) s$use_regularity)

    offset <- F$offset

    internal_bases <- map(conv_smooths, "internal_basis")
    
    lps <- predict(F, type = "lpmatrix")

    regularity_fits <- vector("list", length(conv_smooths))
    
    for (iter in 1:indir.max.iter) {

        if(iter != 1) {
            last_fit <- regularity_fits[[i]]
        }
        offset <- offset * 0
        for(i in 1:length(internal_bases)) {

            if (iter == 1) {
                conv_smooths[[i]]$indirect_fitting_complete <- FALSE
            }
            
            if(conv_smooths[[i]]$indirect_fitting_complete) next
            
            coefs <- coef(F)
            coef_name <- names(coefs)
            ## add a zero intercept
            start_coef_name <- stringr::str_detect(coef_name, conv_smooths[[i]]$term)
            start_coef <- coefs[start_coef_name]

            regularity_model <- conv_smooths[[i]]$regularity_model
            
            active <- internal_bases[[i]]

            if (regularity_model == "distance_gaussian")
            {
                regularity_data <- data.frame(distances = seq(from = 0.01, to = sqrt(2), length.out = 500))
                ob <- s(distances, bs = "bs2", k = active$bs.dim)
                
                Xp <- smoothCon(ob,
                                data = regularity_data,
                                knots = list(distances = active$knots),
                                absorb.cons = TRUE)[[1]]$X
                
                current_term <- c(Xp %*% start_coef)
                ## rescale and recenter regularity model
                rescale <- min(current_term)
                current_term <- current_term - rescale
                regularity_data$current_term <- current_term

                
                gauss_fit <- nls(current_term ~    alpha * exp(- (distances - off )^2 / sigma),
                                 data = regularity_data,
                                 start = list(alpha = max(current_term), sigma = 1, off = 0),
                                 lower = c(0, 0, 0),
                                 algorithm = "port")

                gauss_proj_coef <- lm.fit(Xp, fitted(gauss_fit))$coef

                Xloc <- lps[,coef_name[start_coef_name]]
                noff <- Xloc %*% gauss_proj_coef + rescale 
                
                offset <- offset + c(noff)
            }  


            if(regularity_model == "2d_gaussian") {
                x <- seq(from = -1, to = 1, length.out = 100)
                y <- seq(from = -1, to = 1, length.out = 100)
                regularity_data <- expand_grid(x = x, y = y)

                ob <- te(x, y, k = c(active$margin[[1]]$bs.dim, active$margin[[2]]$bs.dim))
                Xp <- smoothCon(ob,
                                data = regularity_data,
                                absorb.cons = TRUE)[[1]]$X


                current_term <- c(Xp %*% start_coef)
                ## rescale and recenter regularity model
                rescale <- min(current_term)
                current_term <- current_term - rescale
                regularity_data$current_term <- current_term


                gauss_fit <- nls(current_term ~    alpha * exp(- (a * x^2 + b * y^2 + c * x * y) ),
                                 data = regularity_data,
                                 start = list(alpha = max(current_term), a = 1, b = 1, c = 0),
                                 lower = c(1e-5, 1e-5, 1e-5, -10),
                                 upper = c(1e5, 1e5, 1e5, 10),
                                 algorithm = "port")
                

                gauss_proj_coef <- lm.fit(Xp, fitted(gauss_fit))$coef

                Xloc <- lps[,coef_name[start_coef_name]]
                noff <- Xloc %*% gauss_proj_coef + rescale 
                
                offset <- offset + c(noff)

            }
            if (iter != 1) {
                delta <- max(abs(coef(last_fit) - coef(gauss_fit)))
                if (delta < delta.tol) regularity_fits[[i]]$indirect_fitting_complete <- TRUE
            }
            
            regularity_fits[[i]] <- gauss_fit
        }

        F <- updateF(F, G, offset)
        F$regularity_fits <- regularity_fits
    }
    return(F)
}







indirect_bam <- function (
                          formula,
                          data,
                          family,
                          indir.max.iter = 10,
                          delta.tol = 1e-5,

                          ...)
{

    ## G <- gam(formula = formula, family = family, data = data, offset = rep(0, nrow(data)), fit = FALSE)
    G <- bam(formula = formula, family = family, data = data, offset = rep(0, nrow(data)), fit = FALSE, ...)
    F <- bam(G = G)
    class(F) <- c("indirect_gam", class(F))

    smooths <- F$smooth
    conv_smooths <- keep(smooths, \(s)inherits(s, "Convspline.smooth") )
    conv_smooths <- keep(conv_smooths, \(s) s$use_regularity)

    offset <- F$offset

    internal_bases <- map(conv_smooths, "internal_basis")
    
    lps <- predict(F, type = "lpmatrix")

    for (iter in 1:indir.max.iter) {
        if(iter != 1) {
            last_fit <- gauss_fit
        }
        offset <- offset * 0
        for(i in 1:length(internal_bases)) {
            coefs <- coef(F)
            coef_name <- names(coefs)
            ## add a zero intercept
            start_coef_name <- stringr::str_detect(coef_name, conv_smooths[[i]]$term)
            start_coef <- coefs[start_coef_name]

            regularity_model <- conv_smooths[[i]]$regularity_model
            
            active <- internal_bases[[i]]

            if (regularity_model == "distance_gaussian")
            {
                regularity_data <- data.frame(distances = seq(from = 0.01, to = sqrt(2), length.out = 500))
                ob <- s(distances, bs = "bs2", k = active$bs.dim)
                
                Xp <- smoothCon(ob,
                                data = regularity_data,
                                knots = list(distances = active$knots),
                                absorb.cons = TRUE)[[1]]$X
                
                current_term <- c(Xp %*% start_coef)
                ## rescale and recenter regularity model
                rescale <- min(current_term)
                current_term <- current_term - rescale
                regularity_data$current_term <- current_term

                
                gauss_fit <- nls(current_term ~    alpha * exp(- (distances - off )^2 / sigma),
                                 data = regularity_data,
                                 start = list(alpha = max(current_term), sigma = 1, off = 0),
                                 lower = c(0, 0, 0),
                                 algorithm = "port")

                gauss_proj_coef <- lm.fit(Xp, fitted(gauss_fit))$coef

                Xloc <- lps[,coef_name[start_coef_name]]
                noff <- Xloc %*% gauss_proj_coef + rescale 
                
                offset <- offset + c(noff)
            }  


            if(regularity_model == "2d_gaussian") {
                x <- seq(from = -1, to = 1, length.out = 100)
                y <- seq(from = -1, to = 1, length.out = 100)
                regularity_data <- expand_grid(x = x, y = y)

                ob <- te(x, y, k = c(active$margin[[1]]$bs.dim, active$margin[[2]]$bs.dim))
                Xp <- smoothCon(ob,
                                data = regularity_data,
                                absorb.cons = TRUE)[[1]]$X


                current_term <- c(Xp %*% start_coef)
                ## rescale and recenter regularity model
                rescale <- min(current_term)
                current_term <- current_term - rescale
                regularity_data$current_term <- current_term


                gauss_fit <- nls(current_term ~    alpha * exp(- (a * x^2 + b * y^2 + c * x * y) ),
                                 data = regularity_data,
                                 start = list(alpha = max(current_term), a = 1, b = 1, c = 0),
                                 lower = c(1e-5, 1e-5, 1e-5, -10),
                                 upper = c(1e5, 1e5, 1e5, 10),
                                 algorithm = "port")
                

                gauss_proj_coef <- lm.fit(Xp, fitted(gauss_fit))$coef

                Xloc <- lps[,coef_name[start_coef_name]]
                noff <- Xloc %*% gauss_proj_coef + rescale 
                
                offset <- offset + c(noff)

            }


            F$regularity_model <- gauss_fit

            

        }

        F <- updateF(F, G, offset)
    }
    return(F)
}









