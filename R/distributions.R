param <- function(loc, ...) {
 UseMethod('param')
}

param.glm <- function(loc, ...){
    structure(list(object = loc,
              has_predict = TRUE),
              class = 'parameter')
}

param.numeric <- function(loc, ...){
    structure(list(object = loc,
                   has_predict = FALSE),
              class = 'parameter')
}

predict.parameter <- function(object, x){
    if(!object$has_predict){
        if(is.vector(x))
            return(rep(object$object, length(x)))
        if(is.matrix(x) | is.data.frame(x))
            return(rep(object$object, nrow(x)))
    }
    predict(object$object, newdata = x)
}


Normal <- function(mean, sd) {

    if(!is(mean, 'parameter')) mean <- param(mean)
    if(!is(sd, 'parameter')) sd <- param(sd)
    
    structure(list(mean = mean,
                   sd = sd),
              class = c('normal_distribution', 'distribution'))
}


density <- function(distribution, ...){
    UseMethod('density')
}

densityfun <- function(distribution) {
    UseMethod('densityfun')
}
densityfun.distribution <- function(distribution) {
    function(A, X, ...) density(distribution, A, X, ...)
}

density.normal_distribution <- function(distribution, A, X, log = FALSE){
    dnorm(A, predict(distribution$mean, X), predict(distribution$sd, X), log = log)
}

