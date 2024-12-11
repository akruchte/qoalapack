#' Coordinate
#' 
#' coord types are used as locational information in qoalapack and for mediating interaction with a variet of other spatial packages
#' @param x a number
#' @param y a number
#' @return a coordinate representation
#' @export
coord <- function(x, y) {
    if (missing(y)){
        
        if (is_coord(x)) return(x)
        
        stopifnot( matrixlike(x) & ncol(x) >= 2)
        
        y <- x[,2]
        x <- x[,1]
    }
        
    
    new_rcrd(list(x = x, y = y), class = 'coord')
}

matrixlike <- function(ob){
    is.data.frame(ob) | is.matrix(ob) | is.array(ob)
}
    


#' label_coord
#'
#' Add a label to a pre-existing coord object and convert to lcoord type
#' coord x Any -> lcoord
#' @export
label_coord <- function(ob, label) {
    stopifnot(inherits(ob, 'coord'))
    
    new_rcrd(list(x = coordx(ob), y = coordy(ob), label = label), class = c('lcoord', 'coord'))
}

#' Labeled coord
#'
#' Labeled coordinates allow additional metadata be assigned to a coordinate.
#' These inherit generally from coord but allow additional extensions. 
#' @param x a number
#' @param y a number
#' @param label a list or vector of additional metadata
#' @return a coordinate representation with additional labels
#' @details The label must either be the same length as x or y or of length 1 in which case it is recycled to the length of x
#' @export 
lcoord <- function(x = double(), y = double(), label){
    vec_assert(x, ptype = double())
    vec_assert(y, ptype = double())

    new_rcrd(list(x = x, y = y, label = vec_recycle(label, size = vec_size(x))), class = c('lcoord', 'coord'))
}


#' @exportS3Method
vec_ptype_abbr.coord <- function(x, ...){
    'coord'
}

#' @exportS3Method
vec_ptype_abbr.lcoord <- function(x, ...){
    'coord:labeled'
}


#' @exportS3Method
format.lcoord <- function(ob, ...){
    out <- format.coord(ob)
    out <- paste0(out, '@', field(ob, 'label'))
    out
}


#' @exportS3Method
format.coord <- function(ob, ...){
    x <- signif(field(ob, 'x'), 2)
    y <- signif(field(ob, 'y'), 2)
    out <- paste0('(', x, ',', y, ')')
    out
}

#' is_coord predicate
#'
#' coord -> bool
#' check if object inherits from coord
#' @export
is_coord <- function(ob, ...){
    inherits(ob, 'coord')
}


#' @exportS3Method
as_coord.default <- function(ob, ...){
    cl <- class(ob)[[1]]
    vec_cast(ob, coord())
}

#'  Coordinate X
#'
#' coord -> real
#' @param ob a coord
#' @export
coordx <- function(ob){
    field(ob, 'x')
}

#'  Coordinate Y
#'
#' coord -> real
#' @param ob a coord
#' @export
coordy <- function(ob){
    field(ob, 'y')
}


#' @exportS3Method
Norm.coord <- function(ob) {
    sqrt(Norm2(ob))
}

#' @exportS3Method
Norm2.coord <- function(ob) {
    coordx(ob)^2 + coordy(ob)^2
}


#' @export
as_matrix <- function(ob, ...){
    UseMethod("as_matrix")
}
#' @exportS3Method
as_matrix.coord <- function(ob) {
    cbind(coordx(ob), coordy(ob))
}
