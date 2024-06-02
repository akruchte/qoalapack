## generative model candidate library set
## it would be very useful to generalize the model to an entire set of candidate models, retain information about which ones are consistent with the simulation scenario and then prune incompatible models.

## And it would also be nice to be able to set a fixed database where we can use a bayesian bootstrap to approximate the desired model by sampling from that database rather than necessarily having to regenerate everything every single time, at least within a desired level of aproximation accuracy.

## assumptions
positivity: P(A) > 0, P(!A) > 0
estimability: A _|_ X | X, f(U) where f ~ Class
exchangeability: A _|_ Y | X
consistency: Y(A = 1) = Y(1)


## a simulation is an entire data analysis that CAN BE REPEATED arbitrarily many times
simulation <- function( ) {

}


experiment <- function() {}

estimator <- function() {}

estimand <- function() {}

distribution <- function(){}
