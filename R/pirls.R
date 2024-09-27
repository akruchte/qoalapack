x <- rnorm(100)
y <- rpois(100, exp(2 * x + 1) )


irls <- function(y, x, family, w = NULL) {
    
    tol <- 1e-5
    X <- cbind(1, x)
    beta <- lm.wfit(X,y, c(1/(y + 1)))$coefficients
    


    for (i in 1:100) {
        eta <- X %*% beta
        mu <- family$linkinv(eta)
        Z <- eta + (y - mu) /mu


        w <- mu
        beta_new <- lm.wfit(X, Z, w = c(w))$coefficients
        delta <- sum((beta_new - beta)^2)
        beta <- beta_new
        ## if (delta < tol) return(beta)
    }
    return(beta)
    
 }



poisson_grad <- function(beta) {
    mean(x * (y - exp(x * beta)))
}


