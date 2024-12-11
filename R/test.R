if(FALSE){
library(tidyverse)
library(mgcv)


G <- gam(hwy ~ s(cty), data=mpg, fit = FALSE)

ind <- 1:100

G$X <- G$X[ind,]
G$y <- G$y[ind]
G$n <- 100
G$mf <- G$mf[ind,]
G$w <- G$w[ind]
G$offset <- G$offset[ind]
f <- gam(G = G)



f <- gam(G = G, mustart = f$coef)
f2 <- gam(G = G, mustart = f$coef)

G <- gam(hwy ~ s(cty), data=mpg, fit = FALSE)
f3 <- gam(G = G, start = f2$coef)




## example for lasso test



modmat <- matrix(rnorm(1000 * 250), nrow = 250)

funs <- c((function(x) log(x + 100)), sin, cos, tan, abs)

coefs <- rnorm(1000) * rbinom(1000, 1, 0.05)

funs <- sample(funs, size = 1000, replace = TRUE)


Mpred <- modmat
for (i in 1:ncol(modmat)){
    Mpred[,i] <- funs[[i]](modmat[,i])
}


Y <- Mpred %*% coefs

modmat <- as_tibble(as.data.frame(modmat))

data <- bind_cols(modmat, y = Y)

form <- as.formula(paste0('y ~ ', paste0(lapply(1:50, \(ob) paste0('s(V', ob, ')')), collapse = '+')))

fit <- bam(form, data = data, select = TRUE, cluster = nc)
}
