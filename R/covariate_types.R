

#' @export 
evokeable <- function(data) {
    new_vctr(data, 'evokeable')
}
## source('coord.R')

#' @export
evaluate <- function(object, ...){
 UseMethod('evaluate')
}
#' Covariate Placeholder
#'
#' A covariate placeholder type is a vector of coordinates (or other 'lookup' type)
#' which are used to index some arbitrary metadata on a continuous space.
#' In many cases in a spatial setting each point in a data-frame is naturally associated with a local neighborhood
#' of more general information about the spatial problems of interest. However if this information is not something
#' we have necessarily predetermined how to use at each point in advance there is not an obvious way how to represent the data
#' in common data structures.
#'
#' Covariate placeholders store the lookup information directly in a standard vector and use the attributes system to store additional
#' data for runtime/analysis-time lookup. This information is stored locally and used locally.
#'
#' @param coords A vector of coordinates
#' @param data Additional spatial data indexed by coords
#'  @export
covariate_placeholder <- function(data, coords,  ...) {

    if (missing(coords)) {
        stop("default coords not implemented")
        message("Initializing placeholder covariate at default values\n")
        coords <- default_points(extent(data))
    }
    coords <- coord(coords)
    
    pdata <- prepare_placeholder_data(data, ...)
    
    structure(
        coord(coords),
        class = c('covariate_placeholder', class(coords)),
        meta = list(),
        data = pdata$data,
        data_type = pdata$data_type,
        mapping = rep(1, length(coords)))
}

#' @export
prepare_placeholder_data <- function(data, ...){
    UseMethod('prepare_placeholder_data')
}

## covariate placeholder is a coordinate marked with additional information, a context for evaluation, and an evaluation strategy


prepare_placeholder_data.Pcov <- function(object, ...){
    list(data = object, data_type = 'Pcov')
}

#' @exportS3Method
prepare_placeholder_data.ppp <- function(object, ...){

    stop("Not yet implemented for spatstat point process.")

    message("Coercing spatstat point process to Pcov with default arguments.\n" )
    list(data = object, data_type = 'unmarked ppp')
}

#' @exportS3Method
prepare_placeholder_data.im <- function(object, ...){
    list(data = object, data_type = 'image')
}

## Assume all point processes in data will be evaluated at the same set of coords
     ## single   



#' @exportS3Method
c.covariate_placeholder <- function(x, y, ...) {
    if (missing(y)) return(x)
    
    stopifnot(class(y) == 'covariate_placeholder')

    dataeq <- identical(extract_data(x), extract_data(y))

    if (dataeq){
        covariate_placeholder(c(extract_coords(x), extract_coords(y), extract_data(x))) 
    }
    else
    {
    }        
}



#' @exportS3Method
print.covariate_placeholder <- function(object){
    cat('A Covariate Placeholder\n')
    cat('With ',  length(object), ' coordinates.')
    cat('\n')
    NextMethod()
}

## TODO
#' @export
`[.covariate_placeholder` <- function(object, ...){
    warning("NOT IMPLEMENTED FOR METADATA")
    NextMethod()
}

#' extract_data
#' 
#' covariate_placeholder -> metadata
#' @export
extract_data <- function(object) {
  attr(object, 'data')
}

#' extract_data
#' 
#' covariate_placeholder -> coords

#' @export
remap <- function(placeholder, new_coords) {
  covariate_placeholder(extract_data(placeholder), new_coords)
}



