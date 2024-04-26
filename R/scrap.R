source('outcome_models.R')

x <- matrify_ppp(swedishpines, swedishpines, c(128, 128))

distances <- c(fft2shift(x$distances))

orx <- c(x$orx)
ory <- c(x$ory)

angle <- seq(from = -1, to = 1,  length.out = 128 * 2)

angle <- outer(angle, angle, atan2)

angle <- c(angle)

dat <- data.frame(angle = angle, distances = distances)

smooths <- smooth.construct(te(distances, angle,  bs = c('tp', 'cc'), k = c(10, 8)), data = dat, knots = NULL)


basis <- smooths$X

plot_basis <- function(ind) image(matrix(basis[,ind], nrow = 256))


plot_basis(8)


im <-as.matrix(as.im(swedishpines))
buf <- matrix(0, nrow = 256, ncol = 256)
buf[1:128, 1:128] <- im

bmat <- fft2shift(matrix(basis[,1], nrow = 256))


image((Re(fft(fft(bmat) * fft(buf), inverse = TRUE)) / 256^2)[1:128, 1:128])
