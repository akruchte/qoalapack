#' @export 
evokeable <- function(data) {
    new_vctr(data, 'evokeable')
}
## source('coord.R')

#' covariate placeholders should carry the relevant information regarding
#' the appropriate single level entity  information needed in model fitting (such as mgcv::gam)
#'  cases include an age, sex, geo-coordinate (x,y), or possibly higher order coordinates, (x,y,t, w) for extra w
#' @export
placeholder_value <- function(value) {
  structure(value, class = c('placeholder_value', 'numeric'))
}
#' @export
print.placeholder_value <- function(value) {
  cat(str_glue('(({value}))\n\n'))
}

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
covariate_placeholder <- function(data, coords) {
  structure(rep(placeholder_value(1), length(coords)),
            class = c('covariate_placeholder', 'numeric'),
            coords = coords,
            data = data)
}



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
}

## TODO
#' @export
`[.covariate_placeholder` <- function(object, ...){
  covariate_placeholder(data = extract_data(object), coords = extract_coords(object)[...])
}

#' @export
extract_data <- function(object) {
  attr(object, 'data')
}

#' @export
extract_coords <- function(object) {
  attr(object, 'coords')
}


#' extract_data
#' 
#' covariate_placeholder -> coords

#' @export
remap <- function(placeholder, new_coords) {
  covariate_placeholder(extract_data(placeholder), new_coords)
}



