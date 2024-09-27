
param <- function(...) {
    list2(...)
}

gauss_kern <- function(d, param) {
    inter <- param[[1]]
    a <- param[[2]]
    s <- param[[3]]
    
    inter + a * exp(-s * d^2)
}

x <- seq(from = 0, to = 5, length.out = 1000)

smooth <- smooth.construct(s(x, bs = 'bs', k = 7), data = list(x = x), knots = NULL)


find_gaussian <- function(x, y){
    loss <- function(par) mean((gauss_kern(x, par) - y)^2)
    optim(c(1,1,1), loss)$par
}

gk <- gauss_kern(s, list(inter = 1, a = 1, s = 1))


S <- smooth$S
mod <- lm(gk ~ smooth$X)
mod_coef <- coef(mod)[-1]
mod_coef[is.na(mod_coef)] <- 0

gauss_par <- find_gaussian(basis, predict(mod))

best_gaussian <- gauss_kern(s, gauss_par)
fmod <- lm(best_gaussian ~ smooth$X)

gauss_coef <- coef(fmod)[-1]
gauss_coef[is.na(gauss_coef)] <- 0



gp <- gauss_like_pen_coef <- mod_coef - gauss_coef
t(gp) %*% S[[1]] %*% gp

plot(s, basis %*% coef(mod)[-1] + coef(mod)[1])
