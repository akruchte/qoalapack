
get_sampler <- function (dist) {
    dist$sampler
}

Normal <- function(mu = 0, sigma = 1){
    dist <- list(sampler = function ( ) rnorm(1, mu, sigma))
}



Declare <- function(sym, in_env = .GlobalEnv) {
    sym <- as.character(substitute(sym))

    sampler <- function ( ) stop(paste0("Variable ", sym, " called without a binding."))
    variable_evaluation <- function (v) {
        if (missing(v)) return(sampler( ))
        else if (is.numeric(v)) {
            sampler <<- function ( ) v
            return()
        } else {
            sampler <<- get_sampler(v)
        }
    }
    
    makeActiveBinding(sym,
                      variable_evaluation,
                      in_env)

}



Param <- Declare

x <- structure(2, class = "test")
    

    
    declaration_proxy$initial_assignment(distribution)
}

Declare(X)
X <- Normal()




## simple implementation of streams in R

stream_c <- function (x, y) {
    o <- c(x, function ( ) eval(substitute(y)))
    class(o) <- "stream"
    o
}

sc <- function(x,y) stream_c(x,y)

scar <- function(s) {
    s[[1]]
}

scdr <- function(s) {
    s[[2]]()
}

;;stream x num -> stream
`+.stream` <- function (s, y ) {
    force(s)
    stream_c(scar(s) + y, scdr(s)  + y)
}

Nat <- stream_c(1, Nat + 1)

smap <- function(s, f)  {
    sc(f(scar(s)), smap(scdr(s), f))
}


even <- function (x) {
    x %% 2 == 0
}


replicates_of <- function(expr) {
    sc(expr, replicates_of(expr))
}



take <- function(s, n) {
    res <- vector("double", n)
    for (i in 1:n){
        res[i] <- scar(s)
        s <- scdr(s)
    }
    res
}
