if (FALSE) {
library(tidyverse)
x <- seq(from = -1, to = 1, length.out = 250)

y <- x

tbl <- expand_grid(x, y)
r <- seq(from = 0, to = 2, length.out = 250)
 theta <- seq(from = 0, to = 2 * pi, length.out = 360)


tbl <- mutate(tbl, theta = atan2(y, x), r = sqrt(x^2 + y^2))

knots <- seq(from = -1, to = 1, length.out = 25)
knotsR <- seq(from = 0, to = sqrt(2), length.out = 25)
knotsT <- seq(from = -pi, to =  pi, length.out = 25)

desR <- splines::spline.des(knotsR, tbl$r, outer.ok = TRUE)$design
desT <- splines::spline.des(knotsT, tbl$theta, outer.ok = TRUE)$design

comps <- vector("list", 10)

for(i in 1:ncol(desX)) {
    comps[[i]] <- desX * desY[,i]
}

mat <- do.call(cbind, comps)

for(i in 1:ncol(desR)) {
    comps[[i]] <- desT * desR[,i]
}

matRT <- do.call(cbind, comps)



desX <- splines::spline.des(knots, tbl$x, outer.ok = TRUE)$design
desY <- splines::spline.des(knots, tbl$y, outer.ok = TRUE)$design
MpenX <- splines::spline.des(knots, tbl$x, outer.ok = TRUE, derivs = 2)$design
MpenY <- splines::spline.des(knots, tbl$y, outer.ok = TRUE, derivs = 2)$design

dMpenX <- splines::spline.des(knots, tbl$x, outer.ok = TRUE, derivs = 1)$design
dMpenY <- splines::spline.des(knots, tbl$y, outer.ok = TRUE, derivs = 1)$design

d <- sqrt(1/500)

S1 <- (d * t(MpenX) %*% MpenX * d) / 10000
S2 <- (d * t(MpenY) %*% MpenY * d) / 10000
S3 <- d * t(dMpenX) %*% dMpenY * d 

comps <- vector("list", 10)

for(i in 1:ncol(desX)) {
    comps[[i]] <- desX * desY[,i]
}

mat <- do.call(cbind, comps)

coef <- rnorm(ncol(mat))
coefRT <- rnorm(ncol(matRT))

z <- mat %*% coef
z <- matRT %*% coefRT


Z <- matrix(z, nrow = 250)

gaussZ <- with(tbl, exp(-x^2 - y^2))

beta <- c(lm.fit(matRT, gaussZ)$coefficients)

omit <- is.na(beta)
predZ <- matRT[,!omit] %*% beta[!omit]
image(matrix(predZ, nrow = 250), useRaster = TRUE)

}
