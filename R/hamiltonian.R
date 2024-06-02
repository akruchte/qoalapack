library(tidyverse)
library(rlang)
library(vctrs)

dens <- function(x) exp(-(x^2))
metropolis <- function(dens, niter = 100){

    past <- rnorm(1)
    proposal <- rnorm(1)

    samples <- vector('double', niter)

    for (i in 1:niter) {
        increment <- rnorm(1, 0, 1)
        proposal <- past + increment
         A <- pmin(1, dens(proposal) / dens(past))

        proposal <- if_else(runif(1) <= A, proposal, past)
        samples[i] <- proposal
        past <- proposal
    }
    samples
}




if (FALSE) {
    library(tidyverse)
    library(spatstat)
    library(sf)

    library(arrow)


    points <-runifpoint(100)

    write_parquet(coords(points), "~Desktop/test_points.parquet")
}
