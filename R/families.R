## families overwrites default glm and gam families by adding additional functions for calculation of the propensity scores
## other than adding in propensity score calculations, this does not change anything

gaussian <- function(...){
    family <- stats::gaussian(...)
    family$pi <- function(a, mu, sigma){
        dnorm(a, mu, sigma)
    }

    family$mu <- function(object, newdata){
        if(missing(newdata)){
            mu <- fitted(object)
        }
        else {
            mu <- predict(object, newdata, type = 'response')
        }
    }

    family$sigma <- function(object, newdata){
        sd(residuals(fit))
    }
    family
}


scat <- function(...){

}

binomial <- function(...){
    family <- stats::binomial(...)

    family$pi <- function(a, p){
        a * p + (1 - a) * (1 - p)
    }
    family$p <- function(object, newdata){
        predict(object, newdata = newdata, type = 'response')
    }
}
