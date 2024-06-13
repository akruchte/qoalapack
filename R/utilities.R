fft2shift <- function(x) apply(apply(x, 1, pracma::fftshift), 1, pracma::fftshift)
sq <- function(x) x^2
