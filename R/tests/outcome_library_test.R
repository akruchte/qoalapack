library(spatstat)
library(mgcv)
swedishpines
ss <- simulate(ppm(swedishpines))[[1]]

## pass in covariate as image
## basis as kernel, then fft

im <- as.im(swedishpines)
im <- as.matrix(im)
dims <- dim(im)

Q <- quadscheme(swedishpines)
P <- union.quad(Q)
X <- Q$data



Y <- swedishpines
data <- mpl_prepare(Y, Q, list(swedishpines), covariates = NULL, bs = 'bs')

object <- s(x, y, by = ...5)

