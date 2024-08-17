library(spatstat)
library(tidyverse)
library(mgcv)
object <- s(swedishpines, bs = c('conv', 'bs'))
data <- list(swedishpines = swedishpines)
knots <- NULL
