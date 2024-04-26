source('outcome_models.R')


s <- swedishpines


mod <- ppm(s ~ x + y)

mod2 <- ppm(s ~ x + y, interaction = Strauss(0.2))
evalInteraction(s, s, interaction = AreaInter(0.1), correction = 'none')

locs <- coords(s)




bump <- function(d) if_else(d < 1,  exp(-1/(1 - d^2)), 0)

## bump shift initially only supported for circular regions
bump_shift <- function(center, total_intensity, radius){

}

## increase (or decrease) the intensity in a region by some uniform intensity function.
uniform_shift <- function(region, total_intensity){
    
}

## shift intervention is for a continous treatment surface
shift_intervention <- function(){}
regional_intervention <- function(){}
stochastic_intervention <- function(){}
pointwise_intervention <- function(){}


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

more_points_distribution <- function(distribution, npoints, region, inner = TRUE) {
    structure(list(distribution = distribution, npoints = npoints),
              class = c('more_points', 'counterfactual_distribution'))
}

intensity_shift <- function(distribution, shift) {
    structure(list(distribution = distribution, shift = shift),
              class = c('shift', 'counterfactual_distribution'))
}

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
propensity_score.more_points <- function(counterfactual, reference, w, log = TRUE) {
    aw <- area(w)

    propensity_observed <- propensity_score(distribution, log = TRUE)
    propensity_intervention <- -log(aw) * counterfactual$npoints + propensity_score

    propensity_score = -log(aw) * counterfactual$npoints

}

counterfactual_distribution <- function(distribution, counterfactual) {

    structure(list(distribution = distribution,
                   counterfactual = counterfactual),
              class = 'counterfactual_distribution')
}

predict.counterfactual_distribution <- function(object, ...) {
    dist <- object$distribution
    predicted <- predict(dist, ..., type = 'intensity')
    object$counterfactual(predicted)
}

