
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

is_coord <- function(ob, ...){
    inherits(ob, 'coord')
}

as_coord <- function(ob, ...){
    UseMethod('as_coord')
}

as_coord.default <- function(ob, ...){
    cl <- class(ob)[[1]]
    vec_cast(ob, coord())
}

coordx <- function(ob){
    field(ob, 'x')
}

coordy <- function(ob){
    field(ob, 'y')
}


## as_coord.sf


