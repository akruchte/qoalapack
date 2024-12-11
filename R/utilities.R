#' @export
fft2shift <- function(x) apply(apply(x, 1, pracma::fftshift), 1, pracma::fftshift)
#' @export
sq <- function(x) x^2
