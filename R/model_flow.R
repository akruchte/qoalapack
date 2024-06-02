## model semantics
## Let K be a Kahler-manifold object.

## Let T be an arbitrary convex partitioning of space.
## Let Bp be the set of open covers of K such that each element of the cover contains an element of T.
## Now take the set (T \oplus \Bp) \cup (T \ominus Bp)
## where ominus denotes the symmetric difference of the two sets.
## This "topology" like object is naturally understood as being all "components" of the "map" generated
## by graphing the boundaries of all these points.
## From this set we can generate a mesh refinement for calculating likelihoods via the INLA approach.
## However, even recording all the data and fitting the model with the INLA approach may be cost prohibitive. Especially if we want to be able to fit a large number of candidate models and perform model selection, resampling, etc. So we would like to simplify this process further by taking a monte carlo approximation to the inla estimator.

## ;generate swarm of points

## uniform initialization should probably work fine alternative importance initialiazation could put points starting near each observed point, and they should probably be contained by a buffered potential barrier that surrounds the observation space.
## I won't bother making these points interact
## the goal is going to sample points that contribute heavily to the likelihood we are trying to actually estimate. This is precisely the same thing as estimating the model itself.So we will let the hamiltonian approximator be based on the eventual model and then progressively refine the model as we approximate it.
## The goal is to get this process to converge such that at some point we can stop the stochastic dynamic, declare that we will switch to deterministic dynamics, and then progress forwards with expected program behavior.I absolutely have to start using continuations here and racket.





swarm_initialize <- function(){
}

## We want to initialize this by sending the point set out along hamiltonian trajectories that are




## fitting_context -> parset
## fitting_context: (datalike, model, estimator)

## I really do need racket like values here

## I would like to be able to return something like
## return (a,b) where
## mixed semantics
## a <- known-value
## b <- (continuation-expression)
warm_start <- function(fitting_context) {
    ## start the fiting process somehow
    ## initialize the parameter recording structure
}


## we want to pass around model parameters to get restarts, and simplify model fitting without a lot of redundant work


## typical pipeline

## estimand definition
## data_prep
## model_setup
## model_fitting
## model_selection
## model_revision
## model_inference
## estimand estimator
## estimand inference


