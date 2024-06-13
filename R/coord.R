if(interactive()){
    library(vctrs)
}

#' 2d coordinates object
#'
#' @param x A number
#' @param y A number
#'
#' @return coord
#' @export
#'
#' @examples
coord <- function(x = double(), y = double()){
    vec_assert(x, ptype = double())
    vec_assert(y, ptype = double())

    new_rcrd(list(x = x, y = y), class = 'coord')
}


#' Title
#'
#' @param x 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
vec_ptype_abbr.coord <- function(x, ...){
    'coord'
}

format.coord <- function(ob, ...){
    x <- signif(field(ob, 'x'), 2)
    y <- signif(field(ob, 'y'), 2)
    out <- paste0('(', x, ',', y, ')')
    out
}


#' coordinate predicate
#' @export
is_coord <- function(ob, ...){
    inherits(ob, 'coord')
}

#' Convert eligible objects to coords via generic S3 interface.
#' @export
as_coord <- function(ob, ...){
    UseMethod('as_coord')
}

as_coord.default <- function(ob, ...){
    cl <- class(ob)[[1]]
    vec_cast(ob, coord())
}

#'  Coordinate X
#' @export
coordx <- function(ob){
    field(ob, 'x')
}

#' Coordinate Y
#' @export
coordy <- function(ob){
    field(ob, 'y')
}

Norm <- function(ob, ...) {
    UseMethod(ob, ...)
}

Norm2 <- function(ob, ...) {
    UseMethod(ob, ...)
}

Norm.coord <- function(ob) {
    sqrt(Norm2(ob))
}


Norm2.coord <- function(ob) {
    coordx(ob)^2 + coordy(ob)^2
}

## as_coord.sf


