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


## required mgcv function
## this function expects to receive data in the form of a covariate placeholder.
## The levels are included directly in the data
## actual coordinates and pcov data is encoded internally in attributes


#' @export
smooth.construct.conv.smooth.spec <- function(object, data, knots) {

    ## TODO
    ## this should just be a standard Pcov
    conv_data <- purrr::reduce(extract_data(data[[object$term]]), c)
    
    coords <- data[[object$term]]

    extra <- object$xt
    ctxt <- extra$context

    max_dist_prop <- extra$max_dist_prop

    
        ## this should probably be removed
    if (is.null(max_dist_prop)) max_dist_prop <- 0.25

    if (!is.null(ctxt)) {
        initial <- is.na(ctxt$context)
    }
    if (initial) {

        term <- object$term

        ## if no secondary basis is provided default to b-splines
        basis_term <- 'bs'
        if (length(class(object)) > 1) {
            ## TODO restrict to specific eligible options
            basis_term <- str_extract(class(object)[[2]], '[a-zA-Z]+')
        } else {basis_term <- 'bs'}
        

        ## turn into convolution internal basis function
        intern_call <- s(distances,  bs = basis_term, fx = object$fixed, k = object$bs.dim)
        intern_call$label <- paste0('conv(', term, ')')
        
        ## local_data <- list(distances = unique(c(field(conv_data, 'pcov')[[1]]$distances)))
        ## TODO     make user configurable and provide better defaults
        distances <- seq(from = 0, to = max_dist_prop * max(unique(c(field(reduce(conv_data, c), 'pcov')[[1]]$distances))), length.out = 1000)
        neg_buf <- -rev(distances[2:10])
        local_data <- list(distances = c(neg_buf, distances))

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

        return(basis)

    }
}

#' @export
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

    neg_buf <- -rev(distances[2:10])
    local_data <- list(distances = c(neg_buf, distances))


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

    convolve_basis <- function(basis, for_conv, dims, window, coords){
        basis <- fft(basis)
        dim(basis) <- dims * 2


        vec <- Re(fft(basis * for_conv, inverse = TRUE) / prod(dim(basis)))
        spatstat.geom::as.im(vec[1:dims[1], 1:dims[2]], W = window)
    }

    
    for(j in 1:ncols) {

        covar <- lcd$covariate
        dims <- lcd$dims
        window <- lcd$window
        lc <- convolve_basis(Xloc[,j], covar, dims, window, coords)

        interp_basis[[j]] <- lc
       }
        

    ## then in this step apply this to each marked set seperately
    ## in that way everything is now pooled
    basis$interpolation_basis <- list(interp_basis)

    class(basis) <- 'Convspline.smooth'

    basis$X <- Predict.matrix.Convspline.smooth(basis, data)
    basis

}





## adaptive convolutions use a three dimensional convolution
## and then add an additional adaptive surface penalty




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
            interp <- spatstat.geom::interp.im(locinterp[[j]], coords)
            Xmat[row_range, j] <- interp
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




## initialization code taken from mgcv
mgcv_initialize <- function (object, data, knots) {

    ## a B-spline constructor method function
  ## get orders: m[1] is spline order, 3 is cubic. m[2] is order of derivative in penalty.

      if (length(object$p.order)==1) m <- c(object$p.order,max(0,object$p.order-1)) 
      else m <- object$p.order  # m[1] - basis order, m[2] - penalty order


    ## set defaults on m and dim
    if (is.na(m[1])) if (is.na(m[2])) m <- c(3,2) else m[1] <- m[2] + 1
    if (is.na(m[2])) m[2] <- max(0,m[1]-1)
    object$m <- object$p.order <- m
    if (object$bs.dim<0) object$bs.dim <- max(10,m[1]) ## default


    ## knots logic 
    nk <- object$bs.dim - m[1] + 1  # number of interior knots
    if (nk<=0) stop("basis dimension too small for b-spline order")
    if (length(object$term)!=1) stop("Basis only handles 1D smooths")
    x <- data[[object$term]]    # find the data
    k <- knots[[object$term]]
    if (is.null(k))
    {
        xl <- min(x);xu <- max(x)
    }
    else if (length(k)==2)
    { 
        xl <- min(k);xu <- max(k);
        if (xl>min(x)||xu<max(x)) stop("knot range does not include data")
    }
      
      if (!is.null(k)&&length(k)==4&&length(k)<nk+2*m[1]) {
          ## 4 knots supplied: lower prediction limit, lower data limit,
          ##   upper data limit, upper prediction limit
          k <- sort(k)
          dx <- (k[4]-k[1])/(nk-1)
          ko <- c(k[1]-dx*m[1],k[4]+dx*m[1]) ## limits for outer knots
          k <- c(seq(ko[1],k[1],length=m[1]+1),
                 seq(k[2],k[3],length=max(0,nk-2)),
                 seq(k[4],ko[2],length=m[1]+1))
          
      } else if (is.null(k)||length(k)==2) {
          xr <- xu - xl # data limits and range
          xl <- xl-xr*0.001;xu <- xu+xr*0.001;dx <- (xu-xl)/(nk-1) 
          k <- seq(xl-dx*m[1],xu+dx*m[1],length=nk+2*m[1])   
      } else {
          if (length(k)!=nk+2*m[1]) 
              stop(paste("there should be ",nk+2*m[1]," supplied knots"))
      }
      if (is.null(object$deriv)) object$deriv <- 0 
      object$X <- splines::spline.des(k,x,m[1]+1,x*0+object$deriv, outer.ok = TRUE)$design # get model matrix
      if (!is.null(k)) {
          if (sum(colSums(object$X)==0)>0) warning("there is *no* information about some basis coefficients")
      }  
      if (length(unique(x)) < object$bs.dim) warning("basis dimension is larger than number of unique covariates")

}
## modified from mgcv
#' @export

smooth.construct.bs2.smooth.spec <- function(object,data,knots) {
 
  ## now construct derivative based penalty. Order of derivate
  ## is equal to m, which is only a conventional spline in the 
  ## cubic case...
  
  object$knots <- k; 
  class(object) <- "Bspline.smooth"  # Give object a class
  k0 <- k[m[1]+1:nk] ## the interior knots
  object$D <- object$S <- list()
  m2 <- m[2:length(m)] ## penalty orders
  if (length(unique(m2))<length(m2)) stop("multiple penalties of the same order is silly")
  for (i in 1:length(m2)) { ## loop through penalties
    object$deriv <- m2[i] ## derivative order of current penalty
    pord <- m[1]-m2[i] ## order of derivative polynomial 0 is step function
    if (pord<0) stop("requested non-existent derivative in B-spline penalty") 
    h <- diff(k0) ## the difference sequence...
    ## now create the sequence at which to obtain derivatives
    if (pord==0) k1 <- (k0[2:nk]+k0[1:(nk-1)])/2 else {
      h1 <- rep(h/pord,each=pord)
      k1 <- cumsum(c(k0[1],h1)) 
    } 
    dat <- data.frame(k1);names(dat) <- object$term 
    D <- Predict.matrix.Bspline.smooth(object,dat) ## evaluate basis for mth derivative at the k1
    object$deriv <- NULL ## reset or the smooth object will be set to evaluate derivs in prediction! 
    if (pord==0) { ## integrand is just a step function...
      object$D[[i]] <- sqrt(h)*D
    } else { ## integrand is a piecewise polynomial...
      P <- solve(matrix(rep(seq(-1,1,length=pord+1),pord+1)^rep(0:pord,each=pord+1),pord+1,pord+1))
      i1 <- rep(1:(pord+1),pord+1)+rep(1:(pord+1),each=pord+1) ## i + j
      H <- matrix((1+(-1)^(i1-2))/(i1-1),pord+1,pord+1)
      W1 <- t(P)%*%H%*%P
      h <- h/2 ## because we map integration interval to to [-1,1] for maximum stability
      ## Create the non-zero diagonals of the W matrix... 
      ld0 <- rep(sdiag(W1),length(h))*rep(h,each=pord+1)
      i1 <- c(rep(1:pord,length(h)) + rep(0:(length(h)-1) * (pord+1),each=pord),length(ld0))
      ld <- ld0[i1] ## extract elements for leading diagonal
      i0 <- 1:(length(h)-1)*pord+1
      i2 <- 1:(length(h)-1)*(pord+1)
      ld[i0] <- ld[i0] + ld0[i2] ## add on extra parts for overlap
      B <- matrix(0,pord+1,length(ld))
      B[1,] <- ld
      for (k in 1:pord) { ## create the other diagonals...
        diwk <- sdiag(W1,k) ## kth diagonal of W1
        ind <- 1:(length(ld)-k)
        B[k+1,ind] <- (rep(h,each=pord)*rep(c(diwk,rep(0,k-1)),length(h)))[ind]  
      }
      ## ... now B contains the non-zero diagonals of W
      B <- bandchol(B) ## the banded cholesky factor.
      ## Pre-Multiply D by the Cholesky factor...
      D1 <- B[1,]*D
      for (k in 1:pord) {
        ind <- 1:(nrow(D)-k)
        D1[ind,] <- D1[ind,] + B[k+1,ind] * D[ind+k,]
      }
      object$D[[i]] <- D1
    }
    object$S[[i]] <- crossprod(object$D[[i]])
  }
  object$rank <- object$bs.dim-m2  # penalty rank 
  object$null.space.dim <- min(m2)    # dimension of unpenalized space 
 
  object
} ### end of B-spline constructor

#' @export
Predict.matrix.Bspline2.smooth <- function(object,data) {
  object$mono <- 0
  object$m <- object$m - 1 ## for consistency with p-spline defn of m
  Predict.matrix.pspline.smooth(object,data)
}


#' sgam is a thin wrapper around gam that provides revision capabilities. 
## 
#' @export
sgam <- function(formula, data, conv_control, ... ) {
    fit <- gam(formula, data, conv_control, ...)
    class(fit) <- c('sgam', class(fit))
}





## given a certain parametric model use a spline based model to test correctness of the specification
## update the parametric model in one direction or another on the basis of this test

## Once we have a correctly specified test according to the spline approximation we can then revise the spline model to look at max distance specifications or whatnot

spline_spec_test <- function ( ) {
    
}
