source('outcome_models.R')

library(spatstat)

x <- swedishpines

kern <- function(x,y) exp(-5 * sqrt((x^2 + y^2)))

xe <- as.im(x)

ts <- seq(from = -1, to = 1, length.out = 256)

t2 <- seq(from = -2, to = 2, length.out = 256 * 2)
km <- (outer(t2, t2, kern))

zero_pad <- matrix(0, nrow = 256 * 2, ncol = 256 * 2)



f1 <- zero_pad
f1[1:256, 1:256] <- as.matrix(as.im(x, dimyx = 256))
f1 <- fft(f1)

f2 <- fft(apply(apply(km, 1, fftshift), 1, fftshift))


resurf <- fft(f1 * f2, inverse = TRUE) / prod(dim(f1))
resurf <- Re(resurf)[1:256, 1:256]
image(resurf)



y <- c(resurf) + rnorm(length(resurf))

ym <- matrix(y, nrow = 256, 256)
yf <- as.function.im(as.im(ym, window(swedishpines)))
Q <- quadscheme(swedishpines)

dat <- mpl_prepare(swedishpines, quadscheme(swedishpines), ppcov = list(x), dimyx = c(128, 128))

qc <- coords(Q)
yy <- yf(qc[,1], qc[,2])

dat$yy <- yy
fit <- gam(yy ~ s(`...5`, bs = 'conv') - 1, data = dat)

basis <- smoothCon(s(`...5`, bs = 'conv'), data = dat, absorb.cons = TRUE)
basis <- basis[[1]]$X

dat$fitted <- fitted(fit)

smooths <- fit$smooth

smooth_internal <- smooths[[1]]$internal_basis

XX <- PredictMat(smooth_internal, data.frame(distances = seq(from = 0, to = 3, length.out = 400)))

plot(XX[,2:10] %*% coef(fit))



