cos_o <- cos
sin_o <- sin
tan_o <- tan
exp_o <- exp
`+o` <- `+`
`-o` <- `-`
`*o` <- `*`
`^o` <- `^`
`/o` <- `/`

`+.default` <- function(...){
    `+o`
}

`+.dualnum` <- function(...){
    `+o`
}



`-.default` <- function(...){
        `-o`
}

`*.default` <- function(...){
        `*o`
}

`/.default` <- function(...){
        `/o`
}

`^.default`<- function(...){
    `^o`
}

cos <- function(...) {
    UseMethod("cos")
}

cos.dualnum <- function (x) {
    a <- x$primary
    b <- x$dual
    dualnum(cos(a), -sin(a) * b)
}

cos.default <- function (...) {
    cos_o(...)
}

sin <- function(...) {
    UseMethod("sin")
}

sin.dualnum <- function (x) {
    a <- x$primary
    b <- x$dual
    dualnum(sin(a), cos(a) * b)
}

sin.default <- function (...) {
    sin_o(...)
}

tan <- function(...) {
    UseMethod("tan")
}

tan.dualnum <- function (x) {
    a <- x$primary
    b <- x$dual
    dualnum(tan(a), b/(sin(a))^2)
}

tan.default <- function (...) {
    tan_o(...)
}


exp <- function(...) {
    UseMethod("exp")
}

exp.dualnum <- function (x) {
    a <- x$primary
    b <- x$dual
    dualnum(exp(a), b * exp(a))
}

exp.default <- function (...) {
    exp_o(...)
}



`+.dualnum` <- function (x,y, ...) {
    if (is.numeric(x)) {
        return(y + x)
    }
    if (is.numeric(y)){
       return( dualnum(x$primary + y, x$dual))
    }
    if (inherits(y, 'dualnum')) {
        return(dualnum(x$primary + y$primary, x$dual + y$dual))
    }
}

`*.dualnum` <- function (x,y, ...) {
    if (is.numeric(x)) {
        return(y * x)
    }
    if (is.numeric(y)){
        return( dualnum(y * x$primary, y * x$dual))
    }
    if (inherits(y, 'dualnum')) {
        a <- x$primary * y$primary
        b <- x$primary * y$dual
        c <- x$dual * y$primary
        return(dualnum(a, b + c))
    }
}


`/.dualnum` <- function (x,y, ...) {
    if (is.numeric(y)) {
        a <- x$primary
        b <- x$dual
        return(dualnum(a / y, b / y))
    }
    if (is.numeric(x)){
        c <- y$primary
        d <- y$dual
        return( dualnum(x / c, (- x * d) / c^2))
    }
    if (inherits(y, 'dualnum')) {
        a <- x$primary
        b <- x$dual
        c <- y$primary
        d <- y$dual

        return(dualnum(a/c, ((b * c) - (a * d)) / c^2))
    }
}


`^.dualnum` <- function (x,y, ...) {
    if (is.numeric(y)) {
        a <- x$primary
        b <- x$dual
        return(dualnum(a^y, y * a^(y - 1) * b))
    }
    if (is.numeric(x)){
        return(exp(log(x) * y))
    }
    if (inherits(y, 'dualnum')) {
        stop("not yet implemented")
        a <- x$primary
        b <- x$dual
        c <- y$primary
        d <- y$dual

        return(dualnum(a/c, ((b * c) - (a * d)) / c^2))
    }
}



symbolic <- function(symbol) {
    structure(
        list(symbols = substitute(symbol),
             coef = 1L),
        class = 'symbolic'
    )
}

`+.symbolic` <- function(x,y) {
    
}


f <- function(x,y, dind = 0){
    if (any(dind ==1)) {
        x <- dualnum(x)
    }
    if (any(dind ==2)) {
        y <- dualnum(y)
    }
    x * y + 2 * y + 1
}





dualnum <- function (ob, y) {
    if (missing(y)) y <- 1
    structure(list(primary = ob, dual = y), class = 'dualnum')
}
    ## ob <- substitute(ob)
    ## if (is.numeric(ob)) {
    ##     out <- list(primary = ob, dual = y)
    ## }
    ## else
    ## {
    ##     out <- list(primary = ob, dual = y)
    ## }
    ## class(out) <- 'dualnum'
    ## out
## }

eps <- dualnum(0,1)

print.dualnum <- function(ob) {
    a <- ob$primary
    b <- ob$dual

    cat(paste0(a, " + ", b, "ε"))
    cat("\n")
}

`[.dualnum` <- function(dn, num) {
    num
}

dual_ops <- function(f, g) {
    function (x) {
        a <- x[[2]]
        b <- x[[3]]
        structure(list(x[[1]],
                       function(y) f(a(y)),
                       gradient_of(f) (a * b)),
                  class = 'dualnum')
    }
}


## value_of <- function (f) {
##     function (x) {
##         f(x)[[2]]
##     }
## }

## gradient_of <- function (f) {
##     if (typeof(f) == 'closure') {
##         ## call f after "pushing inside" the gradient_call as a symbolic call
##         ## induces autodiff semantics 
##         return (function(...) f(dualnum(...)))
##     }

## if (typeof(f) == 'builtin')
## {
##     if (identical(f, cos)) {
##         return(sin)
##     } else if (identical (f, sin)) {
##         return(function(x) -cos(x))
##     } else if (identical (f, exp)) {
##         return(exp)
##     }
## stop("gradient not yet implemented")
## }
           
## }

    
