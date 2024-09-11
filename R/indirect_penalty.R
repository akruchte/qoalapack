#' Indirect Penalty Calculation
#'
#' Indirectly defined penalties define an additive model penalty in terms of
#' a penalty defined in terms of some (potentially highly non-linear) model.
#' While in principle such a model can be estimated directly the parameters of such a model
#' can be challenging to estimate frequently depending on non-convex losses.
#' They may further require expensive basis set-up at every step of an iterative estimation process.
#' 
#' Further, such models may force excessively strong constraints on the data generating process. 
#' Use of an indirect penalty allows us to specify our models in terms of a practical linear basis expansion which is
#' initialized just once, and penalize the model so as to force it "in the direction of" the underlying motivating model.
#' In this way the correspondence with the motivating model is only approximate. Data driven considerations allow us to fit
#' models that are approximately the same as the motivating model, but encompass much larger classes of functions. 
#'
#' The degree to which the model fits to the regularization model is further defined in terms of itself.
#' We calculate a penalty term (usually a quadratic integral) with respect to the basis used for the model.
#' Then given a set of basis coefficients we find the parameters for the non-linear model which most closely correspond to the
#' parameters of the basis.
#' This model is then taken as the "true" regularity model and the projection of this model onto the basis space is taken to find the best approximation
#' to the regularity model given the basis. This allows the penalty of the regularity model to be calculated approximately and subtracted away from the
#' penalty in the underlying model. 
#' @name indir
#' @export
indirect_penalty <- function(basis, coef, penalty, regularity_model, regularity_control = regularity_model$regularity_control) {

    x <- regularity_control$eval_sequence
    X <- predict(basis, new = x)
    ypred <- X %*% coef

    reg_fit <- regularity_model$fit(x, ypred)
    yreg <- reg_fit$fitted

    reg_coefs <- lm.fit(X, yreg)$coefficients 

    penalty <- evaluate_penalty(coef - reg_coefs, penalty)
}

evaluate_penalty <- function (coef, S = diag(length(coef))) {
    pen <- t(coef) %*% S %*% coef
    pen[[1]]
}


test_data.indirect_gaussian <- function( ) {
    n <- 125

    set.seed(17)
    x <<- seq(from = 0, to = 5, length.out = n)
    y <<- 4 * exp(-x^2 * 2)  +3 + rnorm(n)

    library(mgcv)
    basis <<- smooth.construct(s(x, bs = 'bs'), data = list(x = x), knots = NULL)
}

#' @rdname indir
#' @export
gaussian_regularity <- function( ) {
    structure(list(
        fit = function (){},
        check_valid = function ( ) {},
        regularity_control = list()
    ),
    class = "gaussian_regularity")
}

#' @rdname indir
#' @export
gaussian <- function(d, param) {
    a <- param[[2]]
    s <- param[[3]]
    
    a * exp(-s * d^2)
}


#' @rdname indir
#' @export
find_gaussian <- function(x, y){
    loss <- function(par) mean((gauss_kern(x, par) - y)^2)
    optim(c(1,1), loss)$par
}
