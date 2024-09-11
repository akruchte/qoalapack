#' Splines and Other Basis Functions For Convolutional Spline Models
#'
#' Convolutional spline models depend on the choice of a basis representation with respect to distance.
#' The convolutional basis expansion is then calculated by taking a convolution of the exposure representation with the basis functions.
#'
#' B spline bases are a natural choice and currently the only implemented basis representation for the convolutional spline model.
#'
#' @param basis A b-spline basis object as constructed by `bspline_basis`
#' @param knots A set of knot values (including boundary knots) defining the basis expansion. If missing knots will be defined uniformly.
#' @param n_basis The number of basis functions to be used. A default of 10 is used if this parameter is missing. 
#' @param order The polynomial order of the basis. Convention is to include the intercept in the order so that a cubic polynomial has order 4.
#' @param x A set of points at which the basis expansion is evaluated.
#' @name spline

#' @rdname spline
#' @export
wiggliness_penalty <- function (basis, quadrature_control = "default") {
    quadrature_control <- match.arg(quadrature_control)

    if (quadrature_control == "default") {
        n_quad_points <- 100L
        knots <- unique(sort(basis_knots(basis)))

        n_interval <- length(knots) - 1
        evaluation_points <- vector("list", n_interval)
        deltas <- vector("list", n_interval)
        for (i in 1:n_interval) {
            delta <- (knots[[i + 1L]] - knots[[i]]) / n_quad_points
            evaluation_points[[i]] <- seq(from = knots[[i]], to = knots[[i + 1]], length.out = n_quad_points)
            deltas[[i]] <- rep(delta, n_quad_points)
        }
        
        deriv_points <- do.call(c, evaluation_points)
        integrand_term <- sqrt(do.call(c, deltas))


        deriv_values <- Derivative2(basis, evaluation_points = deriv_points)
        integrand_term <- broadcast(integrand_term, deriv_values)

        integrand <- deriv_values * integrand_term
        S <- t(integrand) %*% integrand
        return(S)
    }    
}

Derivative2 <- function (basis, ..., evaluation_points) {

    if(missing(evaluation_points)){
        evaluation_points <- attr(basis, "x_original")
    }
    knots <- basis_knots(basis)
    order <- basis_order(basis)

    splines::spline.des(knots, evaluation_points, order, derivs = 2, outer.ok = TRUE)$design
}


select_knots <- function (x, n_basis, order)
{

    if (n_basis < 1) {
        n_basis <- 10
    }

    xrange <- range(x)
    
    n_right_boundary_knots <- 2L
    n_left_boundary_knots <- 2L
    
    knots <- seq(from = xrange[1], to = xrange[2], length.out = n_basis)
    interval <- knots[2] - knots[1] 
    knots <- c(-2 * interval, -1 * interval, knots, interval, 2 * interval)

    knots    
}

#' @rdname spline
#' @export
bspline_basis <- function (x,  n_basis = -1, order = 4,
                           knots = select_knots(x = x, n_basis = n_basis, order))
{
    prefab <- splines::spline.des(knots, x, order, derivs = 0, outer.ok = TRUE)

    structure(prefab$design,
              degree = ncol(prefab$design),
              x_original = x,
              knots = prefab$knots,
              order = prefab$order,
              class = "bspline_basis")
}


#' @rdname spline
#' @export
basis_knots <- function(basis) {
    attr(basis, "knots")
}


#' @rdname spline
#' @export
basis_order <- function(basis) {
   attr(basis, "order")
}




broadcast <- function(vec, target_matrix) {
    out <- target_matrix * 0L

    if(! length(vec) == nrow(target_matrix)) {
        stop("Currently only supported for column broadcasting")
    }

    for (col in 1:ncol(out)){
        out[,col] <- vec
    }
    out
}
