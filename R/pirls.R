

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



AdditiveModel <- R6Class("AdditiveModel",
                         public = list(
                             #' Additive Model Constructor
                             #'
                             #' Initializes an additive model
                             initialize = function(basis_like, outcome, family, penalty_model, weights, model_control)
                             {
                                 self$X <- as_basis(basis_like)
                                 self$outcome <- outcome

                                 ## TODO regularity model and regularity parameters
                                 

                                 self$control = model_control
                                 self$weights = weights
                                 
                                 self$beta <- initialize_pirls(self$X, self$outcome)
                             }
                         ),
                         private = list(
                             beta = NULL,
                             basis = NULL,
                             penalty = NULL,
                             penalty_parameter = NULL,
                             control = NULL,
                             weights = NULL
                         )
                         )

                         

initialize_pirls <- function (y, X, S, family, offset, penalty_offset, lambda, w_extern, control = pirls_control()) {
    if (family$family != poisson()) stop("Currently PIRLS only implemented for Poisson family")
     mu <- y + 0.1
    eta <- family$linkfun(mu)

    Z <- eta + (y - mu) /mu
    
    W <- make_weight_matrix(mu * w_extern)

    M1 <- t(X) %*% W %*% X + lambda * S
    M2 <- t(X) %*% W %*% Z

     beta <- solve(M1, M2)
     return(beta)
}

pirls_fit <- function (beta0, y, X, S, offset, subset, penalty_offset, lambda, w_extern, control = pirls_control(), n_step = 1) {

    if (family$family != poisson()) stop("Currently PIRLS only implemented for Poisson family")

    beta <- beta0 + 0.1
    eta <- X %*% beta
    mu <- family$linkinv(eta)
    iter <- 0
    repeat {
        eta <- X %*% beta
        mu <- family$linkinv(eta)
        Z <- eta + (y - mu) /mu

        w <- mu * w_extern
        W <- make_weight_matrix(w)
        
        M1 <- t(X) %*% W %*% X + lambda * S
        M2 <- t(X) %*% W %*% Z

        beta_new <- solve(M1, M2)
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
    return(beta)
}

pirls_step <- function (beta0, y, X, S, subset, offlset, penalty_offset, lambda, w_extern, control = pirls_control(), n_step = 1) {
    if (family$family != poisson()) stop("Currently PIRLS only implemented for Poisson family")

    N <- nrow(X)
    P <- ncol(X)
    
    if (is.null(w_extern)) w_extern <- rep(1, N)

    tol <- control$tol
    max_iter <- control$max_iter

        mu <- y + 0.1
    eta <- family$linkfun(mu)

    Z <- eta + (y - mu) /mu
    
    W <- make_weight_matrix(mu * w_extern)

    M1 <- t(X) %*% W %*% X + lambda * S
    M2 <- t(X) %*% W %*% Z

    beta <- solve(M1, M2)


    for (iter in 1:n_step) {
        eta <- X %*% beta
        mu_old <- mu
        mu <- family$linkinv(eta)
        Z <- eta + (y - mu) /mu

        w <- mu * w_extern
        W <- make_weight_matrix(w)
        XX <- rbind(t(beta), X)
        
        M1 <- t(X) %*% W %*% X + lambda * S
        M2 <- t(X) %*% W %*% Z

        beta <- beta_new
    }
    return(beta)
}

pirls_fit_with_folds <- function ( ) {

}

pirls_fit_with_spatial_folds <- function ( ) {

}





initialize_regularity_model <- function ( ) { }

update_regularity_model <- function ( ) { }



get_offset_from_regularity_model <- function ( ) { }
get_penalty_from_regularity_model <- function ( ) { }
regularity_model_initialize_penalty_parameter <- function ( ) { }


regularity_model_model_projection <- function ( ) { }
fit_regularity_model <- function ( ) { }



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

#' Spatial Covariate
#'
#' spatial covariate
spatial <- function (var)
{

}

spatial_basis <- function ( var )
{
    
}


matern_spde <- function ( var )
{

}
