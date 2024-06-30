library(Rcpp)
library(here)

library(spatstat)

x <- rpoispp(10)


f <- Rcpp::sourceCpp(file = here("R/src/test.cpp"))
