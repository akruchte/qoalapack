library(tidyverse)

library(RandomFields)
library(mgcv)

set.seed(13)


generating_model <- RMmatern(0.9, var = 3) + RMnugget(var = 0.03)

x <- seq(from = -1, to =1, length.out = 128)
field <- as.matrix(RFsimulate(generating_model, x,x))

cov1 <- as.matrix(RFsimulate(RMmatern(3, var = 1), x, x))
cov2 <- as.matrix(RFsimulate(RMmatern(0.1, var = 1), x,x ))

exposure <- 2 * cov1 - cov2 + field 

## ecuts <- quantile(exposure, probs = seq(0, 1, length.out = 7))
## exposure <- cut(exposure, breaks = ecuts, ordered = TRUE)

dat <- expand_grid(x = x, y = x)
dat <- mutate(dat, exposure = c(exposure), cov1 = c(cov1), cov2 = c(cov2))

fit <- gam(exposure ~ cov1, data = dat)


