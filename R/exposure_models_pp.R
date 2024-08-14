bump <- function(d) if_else(d < 1,  exp(-1/(1 - d^2)), 0)

## bump shift initially only supported for circular regions
#' @export
bump_shift <- function(center, total_intensity, radius){

}

## increase (or decrease) the intensity in a region by some uniform intensity function.
#' @export
uniform_shift <- function(region, total_intensity){
    
}

## shift intervention is for a continous treatment surface
#' @export
shift_intervention <- function(){}
#' @export
regional_intervention <- function(){}
#' @export
stochastic_intervention <- function(){}
#' @export
pointwise_intervention <- function(){}


#' @export
intervene <- function(exposure, intervention, ...){
    UseMethod('intervene')
}
intervene.more_points <- function(...){
    
}


## propensity_score <- function(distribution, intervened, reference_measure = function(...) 1, log = TRUE)
## {
##     ## hacky bad implementation right now only for ppp
##     sum(log(predict(distribution, locations = coords(intervened)))) - integral(predict(distribution))
## }
#' @export
more_points_distribution <- function(distribution, npoints, region, inner = TRUE) {
    structure(list(distribution = distribution, npoints = npoints),
              class = c('more_points', 'counterfactual_distribution'))
}

#' @export
intensity_shift <- function(distribution, shift) {
    structure(list(distribution = distribution, shift = shift),
              class = c('shift', 'counterfactual_distribution'))
}


#' @export
propensity_score.shift <- function(counterfactual, reference, observed, log = TRUE){
    ## For an additive shift poisson distribution the normalizing factors cancel each other out in the propensity score and need not be computed, except for the factor of the integral of the shift
    part_norm <- log(integral(counterfactual$shift))

    num <- sum(log(predict(counterfactual, locations = observed, type = 'intensity')))
    denom <- sum(log(predict(counterfactual, locations = observed, type = 'intensiy')))

    pscore <- part_norm + num - denom

    if(log) return(pscore)
    else return(exp(pscore))
}

## for the point process that assigns n additional points independently to a point process
## the total density of the process is the density of the observed poisson process times the density of each
## of the additional uniform point processes
## assuming these are uniform over the observation window, this is simply a multiplicative factor by (1/|w|)^n

#' @export
propensity_score.more_points <- function(counterfactual, reference, w, log = TRUE) {
    aw <- area(w)

    propensity_observed <- propensity_score(distribution, log = TRUE)
    propensity_intervention <- -log(aw) * counterfactual$npoints + propensity_score

    propensity_score = -log(aw) * counterfactual$npoints

}

#' @export
counterfactual_distribution <- function(distribution, counterfactual) {

    structure(list(distribution = distribution,
                   counterfactual = counterfactual),
              class = 'counterfactual_distribution')
}


#' @export
predict.counterfactual_distribution <- function(object, ...) {
    dist <- object$distribution
    predicted <- predict(dist, ..., type = 'intensity')
    object$counterfactual(predicted)
}

