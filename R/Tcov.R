## Tcov should represent covariates that are represented by GIS vector type objects
## in particular, polygons are likely the only type of representation that makes sense in this context


## for example a tesselation

## Tcov is similar in spirit to Rcov, but handles data via a functional representation
## this may be under some circumstances more amenable to handling uncertainty due to measurement error
## Tcov requires a basis or model for the representation of the data
Tcov <- function(..., formula, method, family){
    object <- list2(...)
    fitted <- lapply(object, function(ob) Tcov_prepare(ob, ...))

    new_rcrd(list(fitted = fitted), class = 'Tcov')
}

format.Tcov <- function(ob, ...){
    rep('Tcov', vec_size(ob))
}


## `+.Tcov` <- function(ob, ob2) {
##     if (inherits(ob2, 'numeric')) {
##         result <- map2(field(ob, 'fitted'),
##              ob2,
##              function(f1, f2) function(x,y) f1(x,y) + f2)
##         return(Tcov(!!!result))
##     }
##     if(inherits(ob2, 'Tcov')){

##         result <- map2(field(ob, 'fitted'),
##              field(ob2, 'fitted'),
##              function(f1, f2) function(x,y) f1(x,y) + f2(x,y))
##         return(Tcov(!!!result))
##     }
##     stop('Not implemented')
## }


## `*.Tcov` <- function(ob, ob2) {
##     if (inherits(ob2, 'numeric')) {
##         result <- map2(field(ob, 'fitted'),
##              ob2,
##              function(f1, f2) function(x,y) f1(x,y) * f2)
##         return(Tcov(!!!result))
##     }
##     if(inherits(ob2, 'Tcov')){

##         result <- map2(field(ob, 'fitted'),
##              field(ob2, 'fitted'),
##              function(f1, f2) function(x,y) f1(x,y) * f2(x,y))
##         return(Tcov(!!!result))
##     }
##     stop('Not implemented')
## }

## `-.Tcov` <- function(ob, ob2) {
##     if (inherits(ob2, 'numeric')) {
##         result <- map2(field(ob, 'fitted'),
##              ob2,
##              function(f1, f2) function(x,y) f1(x,y) - f2)
##         return(Tcov(!!!result))
##     }
##     if(inherits(ob2, 'Tcov')){

##         result <- map2(field(ob, 'fitted'),
##              field(ob2, 'fitted'),
##              function(f1, f2) function(x,y) f1(x,y) - f2(x,y))
##         return(Tcov(!!!result))
##     }
##     stop('Not implemented')
## }

## `/.Tcov` <- function(ob, ob2) {
##     if (inherits(ob2, 'numeric')) {
##         result <- map2(field(ob, 'fitted'),
##              ob2,
##              function(f1, f2) function(x,y) f1(x,y) / f2)
##         return(Tcov(!!!result))
##     }
##     if(inherits(ob2, 'Tcov')){

##         result <- map2(field(ob, 'fitted'),
##              field(ob2, 'fitted'),
##              function(f1, f2) function(x,y) f1(x,y) / f2(x,y))
##         return(Tcov(!!!result))
##     }
##     stop('Not implemented')
## }




evaluate.Tcov <- function(ob, locations){
    x <- coordx(locations)
    y <- coordy(locations)
    preds <- map(field(ob, 'fitted'), function(f) f(x,y))
    reduce(preds, c)
}



Tcov_prepare <- function(object, ...){
     UseMethod('Tcov_prepare')
}

Tcov_prepare.function <- function(object, ...) {
    stopifnot(names(formals(object)) == c('x', 'y'))
    object
}

Tcov_prepare.data.frame <- function(data, formula, method, family){
    method(formula, family = family, data = data)
}
