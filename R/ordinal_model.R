## this is horribly incorrect, but kept for its own sake

data(wine, package = 'ordinal')
library(tidyverse)
data_subsets <- function(cat) {
    ls <- levels(cat)
    baselevel <- ls[1]
    others <- ls[-1]

    which.base <- which(cat == baselevel)

    subsets <- vector('list', length = length(others))

    for (i in 1:length(subsets)){
        subsets[[i]] <-which(cat == others[i])
    }
    list(base = which.base, subsets = subsets)

}


intercept_pass <- function(data, subsets, coefficients = rep(1, length(subsets))){

    base <- subsets$base
    subsets <- subsets$subsets
    for (i in 1:length(subsets)) {
        subset <- c(base, do.call(c, subsets[1:i]))
        local_data <- data[subset,]

        fit <- glm.fit(rep(1, nrow(local_data)),
                       local_data[['...loutcome...']],
                       local_data[['beta_offset']],
                       start = max(coefficients[(i-1):i]),
                       family = binomial())
        coefficients[[i]] <- fit$coefficient
    }
    coefficients
}

make_estfun <- function(X, data, subsets) {
    subs <- lapply(subsets$subsets, function(sub) do.call(c, list(subsets$base, sub)))
    Xs <- lapply(subs, function(sub) matrix(X[sub,], ncol = ncol(X)))
    outcomes <- lapply(subs, function(sub) data[sub, '...loutcome...'])
    
    function(beta, intercepts) {
        est <- 0
        for (i in 1:length(subsets)) {
            X <- cbind(1, Xs[[i]])
            Y <- outcomes[[i]]
            estcom <- t(X) %*% (Y - plogis(X %*% c(intercepts[[i]], beta)))
            est <- est + (t(estcom) %*% estcom)[[1]]
        }
        est
    }
}


parameter_pass <- function(estfun, intercepts, beta){
    to_optim <- function(beta) estfun(beta, intercepts)

    optim(beta, to_optim)
}

estimate_model <- function(data, X, outcome_var, max.iter = 100, tol = 1e-8) {

    base_level <- levels(data[[outcome_var]])[[1]]
    data[['...loutcome...']] <- if_else(data[[outcome_var]] == base_level, 0L, 1L)
    subsets <- data_subsets(data[[outcome_var]])
    estfun <- make_estfun(X, data, subsets)

    initialize <- glm.fit(cbind(1, X), data[['...loutcome...']], family = binomial())$coefficient
    intercepts <- rep(initialize[[1]], length(subsets$subsets))
    beta <- initialize[-1]
    data$beta_offset <- X %*% initialize[-1]
    intercepts <- intercept_pass(data, subsets, coefficients = intercepts)
    
    for (i in 1:max.iter){
        old_par <- c(intercepts, beta)
        data$beta_offset <- X %*% beta
        intercept_pass(data, subsets, coefficients = rep(initialize[[1]], length(subsets$subsets)))
        beta <- parameter_pass(estfun, intercepts, beta)$par

        new_par <- c(intercepts, beta)
        if(max(abs(new_par - old_par)) < tol)
            return(list(intercepts = intercepts,
                        beta = beta, iterations = i))
    
    }
    return(list(intercepts = intercepts, beta = beta))

    
}

for_pass <- mutate(wine, ...loutcome... = if_else(rating == levels(rating)[[1]], 0, 1))



