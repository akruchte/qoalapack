library(testthat)
library(disslib)

test_check("disslib")

library(spatstat)
p <- Pcov(swedishpines, W = swedishpines$window, dimyx = c(128, 128))


