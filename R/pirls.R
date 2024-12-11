## main function fitting loop
fit <- function ( ) {
    parameter_context

    ## we will update the global model with all possible information,
    ## then given the global model revise penalty parameters (estimated with CV)
    ## updating regularity models occurs while updating penalties
    ## updates to model parameters happens inside these loops
    
    
    update_global()
    update_penalty()

    update_regularity()

    update_internal()
}



#' Penalized Iterative Least Squares Control Parameters
#'
#' pirls_control creates a control parameter list.
#' Default arguments for all necessary control parameters are provided.
#' This function can be called with keyword arguments to specify alternative control
pirls_control <- function(..., tol = 1e-5, max_iter = 10) {
    structure(list(tol = 1e-5, max_iter = 10, ...),
              class = "pirls_control")
}

#' @exportS3Method
print.pirls_control <- function(object) {
    nms <- names(object)
    for (i in seq_along(object)) {
        fmt_string <- paste0(nms[[i]], ": ", object[[i]], "\n")
        cat(fmt_string)
    }
    return(invisible(object))
}

make_weight_matrix <- function(w) {
    return(diag(c(w)))
}
                         

initialize_pirls <- function (y, X, family, w_extern) {
    if (family$family != "poisson") stop("Currently PIRLS only implemented for Poisson family")
  
     mu <- y + 0.1
    eta <- family$linkfun(mu)

    Z <- eta + (y - mu) /mu
    
    W <- make_weight_matrix(mu * w_extern)
  

    M1 <- t(X) %*% W %*% X + diag(1, ncol(X))
    M2 <- t(X) %*% W %*% Z

  beta <- solve(M1, M2)
     return(beta)
}

#lambda is a vector of penalty parameters
#S is a list of penalty matrices with possibly regularity models
#' @export
make_penalty <- function(S, lambda, n_unpenalized) {
  
  S_construct <- S
  penalty_matrix_size <- 0
  for (mat_index in 1:length(S)) {
     S_construct[[mat_index]] <- S[[mat_index]] * lambda[mat_index]
     penalty_matrix_size <- penalty_matrix_size + ncol(S[[mat_index]])
  }
  penalty_matrix_size <- penalty_matrix_size + n_unpenalized
  Sout <- diag(0, penalty_matrix_size)
  start <- n_unpenalized + 1L
  for (i in 1:length(S)){
    nc <- ncol(S_construct[[i]]) 
    indices <- start:(start + nc - 1)
    start <- start + nc
    Sout[indices, indices] <- S_construct[[i]]
  }
  Sout
}

#' regularity models
#' each regularity model is a list with evaluation points, basis at those evaluation points,
#' model to be fit to, offset prediction, and coefficient mapping


## penalty_construct indicates mapping indices between beta and the penalty matrices.
pirls_fit <- function (usebeta0,
                      betaunpen,
                       betas,
                       y,
                       Xunpen,
                       Xpen, S, 
                       subset, lambda,
                       penalty = TRUE,
                       w_extern = 1,
                       pcontrol = pirls_control(), family = poisson())
{
  
    if (family$family != "poisson") stop("Currently PIRLS only implemented for Poisson family")
    
    X <- cbind(Xunpen, do.call(cbind, Xpen))
    
    #precondition check
      NTEST <- length(y)
      stopifnot(nrow(X) == NTEST)
      stopifnot(length(w_extern) == 1 | length(w_extern) == NTEST)
      
      PTEST <- length(S)
      stopifnot(length(betas) == PTEST)
      stopifnot(length(lambda) == PTEST)
      stopifnot(length(Xpen) == PTEST)
      stopifnot(length(regularity_models) == PTEST)
    
    tol <- pcontrol$tol
    max_iter <- pcontrol$max_iter
    
    if(usebeta0){
      beta <- c(betaunpen, do.call(c, betas))
    } else {
      beta <- initialize_pirls(y, X, family = family, w_extern = w_extern)

      betaunpen <- beta[1:length(betaunpen)]
      start_index <- length(betaunpen) + 1
      for (i in 1:length(betas)) 
      {
        npar_to_grab <- length(betas[[i]])
        betas[[i]] <- c(beta[start_index:(start_index + npar_to_grab - 1)])
        start_index <- start_index + npar_to_grab
      }
    }

    eta <- X %*% beta
    mu <- family$linkinv(eta)

    #pack together penalties for fitting
    if (penalty) {
      lambdaS <- make_penalty( S, lambda, ncol(Xunpen))
    } else {
        lambdaS <- 0      
    }
    
    #primary fitting loop
    iter <- 0
    repeat {


    start_index <- length(betaunpen) + 1

        eta <- X %*% beta 
        mu <- family$linkinv(eta)
        Z <- eta + (y - mu) /mu

        w <- mu * w_extern * subset
        W <- make_weight_matrix(w)
        
        
        M1 <- t(X) %*% W %*% X + lambdaS
        M2 <- t(X) %*% W %*% Z

        beta_new <- solve(M1, M2)     

    for (i in 1:length(betas)) 
    {
      npar_to_grab <- length(betas[[i]])
      betas[[i]] <- c(beta[start_index:(start_index + npar_to_grab - 1)])
      start_index <- start_index + npar_to_grab
    }


    delta_beta <- sum((beta - beta_new)^2)
    beta <- beta_new

    if (delta_beta < tol) {
        stop_condition <- "tolerance met"
         break
    }

    iter <- iter + 1
    if (iter >= max_iter) {
         stop_condition <- "max pirls iterations"
    }

    }
    
    mse <- mean((y[which(subset == 0L)] - family$linkinv(X[which(subset == 0L),] %*% beta))^2)
  
#    unpack beta
    betaunpen <- beta[1:length(betaunpen)]
    return(list(betaunpen = betaunpen, betas = betas, mse = mse))
}

lambda_select <- function(folds, betaunpen,
                          betas,
                          y, Xunpen, Xpen, S, lambda0, 
                          regularity0, 
                          w_extern = 1, control = pirls_control(), family = poisson()) 
{
  loss <- function(LAMBDA){
    
    mses <- pirls_fit_with_folds(folds = folds, beta0 = beta0, y = y, Xunpen = Xunpen, Xpen = Xpen, 
                                S = S, 
                               lambda0 = LAMBDA, w_extern = w_extern, control = control, family = family)$mses
    mean(mses)
  }
    
    optim(lambda0, loss, lower = 1e-5)
}

#' Outer iteration to find lambda values
#' @export
pirls_fit_with_folds <- function (folds,
                                  betaunpen, 
                                  betas,
                                  y, Xunpen,Xpen, S,
                                  regularity,
                                  lambda0,
                                  w_extern = 1,
                                  pcontrol = pirls_control(),
                                  family = poisson()) 
{
  betaunpens <- vector("list", ncol(folds))
  betapens <- vector("list", ncol(folds))
  mses <- vector("double", ncol(folds))

  for (fold in 1:ncol(folds)){
    current_fold <- c(folds[,fold])
    
    active <- pirls_fit(
      usebeta0 = TRUE,
      betaunpen = betaunpen,
                        betas = betas,
                         y = y,
                         Xunpen = Xunpen, 
                         Xpen = Xpen, 
                         S = S,
                         subset = current_fold,
                         lambda = lambda0,
                         w_extern = w_extern,
                         pcontrol = pcontrol,
                         penalty = TRUE,
                         family = family)
    
    betaunpens[[i]] <- active$betaunpen
    betapens[[i]] <- active$betas
    mses[fold] <- active$mse
  }
  
  return(list(betaunpens = betaunpens, betapens = betapens, mses = mses))
}

pirls_fit_with_spatial_folds <- function ( ) {

}


#' Convolutional Basis Representation
#'
#' Setup a convolutional basis (lazily) for use in a convolutional spline model
#' This function is used for setting up model specifications for spatial models with spatial covariates.
conv <- function (var, basis = "bs2",
                  method = c("fourier", "direct"),
                  parametric = FALSE,
                  resolution,
                  max_distance,
                  mark_model,
                  interactions,
                  regularity_model)
{
    
}

#' @export
simple_test <- function () {
  set.seed(100)
  y <- rpois(1000, 1)
  X <- cbind(1, rnorm(1000))
  beta0 <- rnorm(2)
  
  lambda <- list(0)
  S <- list(diag(1, 2))
  penalty_construct <- list(nparams = 2, mapping = list(c(1,2)))
  
  beta0 <- pirls_fit(beta0, y, X, S, penalty_construct, subset = 1, lambda = list(0))
  
  folds <- random_folds(X, 5)
  lambda_select <- lambda_select(folds, beta0, y, X, penalty_construct = penalty_construct, S = S, lambda0 = 0, w_extern = 1)

}


