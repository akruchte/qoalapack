#' @export
matern <- function(d, sig2, v, rho) {
   dplyr::if_else (d == 0,
             1.0,
             sig2 * 2 ^ (1 - v) / gamma(v) * (sqrt(2 * v) * d / rho)^v * besselK(sqrt(2 * v) * d / rho, v))
}

#' @export
sim_matproc <- function() {
    x <- seq(from = -5, to = 5, length.out = 100)
    dist <- function(x, y) sqrt(x^2 + y^2)
    matkern <- function(x,y) matern(dist(x,y), 1, 1, 1)
    distmat <- outer(x, x, matkern)


    noise <- distmat
    noise[] <- rnorm(nrow(noise) * ncol(noise))

    matproc <- Re(fft(fft(noise) * fft(distmat), inverse = TRUE))
    matproc
}
