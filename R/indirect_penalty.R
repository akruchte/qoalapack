## NEXT STEPS


## set up the loop of setting up the automatic selection of penalty parameters fitting the indirect gaussian penalty and so on

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

#' @rdname indir
#' @export
gaussian_kernel <- function(d, param) {
    a <- param[[1]]
    s <- param[[2]]
    
    a * exp(-s * d^2)
}


#' @rdname indir
#' @export
find_gaussian <- function(x, y){
    loss <- function(par) mean((gaussian_kernel(x, par) - y)^2)
    optim(c(1e-5,1e-5), loss)$par
}

#' @export
make_gaussian_regularity_structure <- function(convolutional_basis) {
  distance_basis <- convolutional_basis$distance_basis
  knots <- basis_knots(distance_basis)
  
  ev_points <- seq(from = 0, to = max(knots), length.out = 1000)
  gauss_fit_basis <- bspline_basis(ev_points, n_basis = length(knots), knots = knots)
  list(x = ev_points, 
       basis = gauss_fit_basis,
       par = rep(1e-5, ncol(distance_basis)))
}

#' @export
update_gaussian_regularity <- function(rs, model_coef) {
  update <- fit_gaussian_regularity(rs$x, rs$basis, y = rs$basis %*% model_coef)
  rs$par <- update
  rs
}

#' @export
gaussian_regularity_offset <- function(rs, bs) {
  bs %*% rs$par
}

#' @export
fit_gaussian_regularity <- function(x, basis, y) {
  gaussian_fit <- find_gaussian(x,y)
  gaussiany <- gaussian_kernel(x, gaussian_fit)
  
  regularity_fit <- lm.fit(basis, y)$coefficients
  c(regularity_fit)
}

