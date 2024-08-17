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






## continuous

## numerator model is used for stabilization weights
## it takes in an exposure model as an argument and fits
## the model required for propensity score stabilization
numerator_model <- function(exposure_model){

}
## exposure models need to receive data in the form of a marked process, but where the marks are considered the outcomes rather than the locations
## The basic quadrature scheme tools can then be used, but used for the sake of point evaluations of the covariates at the site observations
## exposure should be represented as a marked process, the marks are assumed to be the response values for modeling.
exposure_model <- function(formula, data, family, method = 'gam', ..., numerator_model =TRUE){ 
    exposure <- data[[as.character(formula[[2]])]]
    
    coords <- as_coord(exposure)
    response_values <- values(exposure)

    covs <- data
    covs[[as.character(formula[[2]])]] <- NULL
    quad <- evaluate_quadrature(covs, at = coords, formula = formula)
    quad[[as.character(formula[[2]])]] <- response_values
    model_engine <- select_engine(method)


    mod_fit <- model_engine(formula, data = data, family = family, ...)

    mod_fit$quad_locations <- coords

    mod_fit
}



select_engine <- function(method) {
    switch(method,
           gam = gam,
           glm = glm)
}


## evaluate each component of a spatial data frame in the required ways to construct the necessary quadrature scheme
evaluate_quadrature <- function(data, at, ..., formula) {
    map_dfc(data, function(ob) evaluate(ob, locations = at))
}

RUN <- FALSE

if(RUN){

    dat <- read_rds('../data_agg/_targets/objects/pm25_cube')

    tdat <- as_tibble(dat) |>
        filter(time %in% time[1])

    tdat <- drop_na(tdat)
    fit <- gam(pm25 ~ s(lon, lat, bs = 'gp'), family = scat(), data = tdat)
    
}


## allow numerator model to be prefit
## if not prefit, use the type of numerator model provided
## if not provided at all, use a sane default
propensity_fun <- function(propensity_model, prediction_locations, numerator_model) {

    outcome_term <- as.character(propensity_model$formula[[2]])
    
    linkinv <- propensity_model$family$linkinv
    family <- propensity_model$family$family
    
    is_binomial <- family %in% c('binomial')
    is_ordinal <- family %in% c('ocat')
    is_continuous <- !is_binomial & !is_ordinal

    mu <- predict(propensity_model, newdata = prediction_locations)
    marginal_outcome <- propensity_model$y

    
    if (is_binomial) {

        prop1 <- mean(marginal_outcome)
        p <- linkinv(mu)
        propfun <- function(A) {
            prop1 / (p * A + (1 - p) * (1 - A))
        }

    }

    if(is_continuous){
        if(missing(numerator_model)) {
            ## use default marginal normal model
            numerator_mu <- mean(marginal_outcome)
            numerator_sd <- sd(marginal_outcome)
            numerator_model <- function(A) dnorm(A, numerator_mu, numerator_sd)
        }

        denominator_mu <- linkinv(mu)
        ## homoskedastic normal case

        denominator_sd <- sd(propensity_model$residuals)
        densfun <- dnorm
        
        propfun <- function(A) {
            numerator_model(A) / dnorm(A, denominator_mu, denominator_sd)
        }
    }

    return(propfun)
}

## fit <- bam(pm25 ~ s(lon, lat, bs = 'ad', k = 150), data = tdat)
## pp <- fitted(fit)
## tdat$pp <- pp

## tdat$pp <- pp  C-c C-c
