x <- seq(from = -1, to = 1, length.out = 128)
distr <- outer(x, x, function(x,y) sqrt(x^2 + y^2))


knots <- log1p(seq(from = 0, to = 1, length.out = 10))

bs <- splines::spline.des(knots, c(distr), outer.ok = TRUE)


p <- function(i) image(matrix(bs$design[,i], nrow = 128))
