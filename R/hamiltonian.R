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
