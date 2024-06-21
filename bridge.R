library(tidyverse)
library(bridgesampling)
library(spatstat)


win <- owin()

randfun <- function() rnorm(1)

dens <- function(x) dnorm(x)

old <- randfun()
proposal <- randfun()

prob <- runif(1)
if (min(dens(proposal) / dens(old), 1) < prob){
    old <- proposal
}



truth_mod <- ppm(proposal ~ x + y)

p <- predict(truth_mod)
ratio <- 
ratio <- log(predict(truth_mod)) - log(predict
## simple sampler

