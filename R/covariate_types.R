## covariate representations should have methods that make them easily coercible for representation in modeling
## specifically, there should be a method for evaluating the covariate representation at the points of a quadrature scheme

#' @export
format.spatial_covariate <- function(object, ...){
    type <- class(object)[[1]]
    rep(type, vctrs::vec_size(object))
}


strip_meta <- c

#' @export 
evokeable <- function(data) {
    new_vctr(data, 'evokeable')
}
## source('coord.R')

#' @export
evaluate <- function(object, ...){
 UseMethod('evaluate')
}


#' @export
evaluate.spatial_covariate <- function(object, ...){
    stop('Evaluate method not implemented')
}
## @param ... <[`dynamic-dots`][rlang::dyn-dots]> What these dots do.


placeholder_recycle <- function(object, target) {

}


#' @export
plot.spatial_covariate <- function(cov, W, ...){
    plot(spatstat.geom::as.im(function(x,y) evaluate(cov, coord(x,y), ...), W))
}

#' @export
as.im.spatial_covariate <- function(cov, W, ...){
        spatstat.geom::as.im(as_Fcov(function(x,y) evaluate(cov, coord(x,y), ...), W))
}

## @param ... <[`dynamic-dots`][rlang::dyn-dots]> What these dots do.


#' covariate placeholders should carry the relevant information regarding
#' the appropriate single level entity  information needed in model fitting (such as mgcv::gam)
#'  cases include an age, sex, geo-coordinate (x,y), or possibly higher order coordinates, (x,y,t, w) for extra w



## the goal of these is to get mgcv to work in a friendly way
## ' min.covariate_placeholder <- function(pl, ...) {min(extract_coords(pl)$x, ...)}
## ' max.covariate_placeholder <- function(pl, ...) {max(extract_coords(pl)$x, ...)}]
#'  @export
covariate_placeholder <- function(coords, data, ...) {
    pdata <- prepare_placeholder_data(data, ...)
    
    structure(
        coords,
        class = c('covariate_placeholder'),
        meta = list(),
        data = pdata$data,
        data_type = pdata$data_type,
        mapping = rep(1, length(coords)))
}

prepare_placeholder_data <- function(data, ...){
    UseMethod('prepare_placeholder_data')
}

## covariate placeholder is a coordinate marked with additional information, a context for evaluation, and an evaluation strategy
#' @exportS3Method
prepare_placeholder_data.ppp <- function(object, ...){
    list(data = object, data_type = 'unmarked ppp')
}

prepare_placeholder_data.ppp <- function(object, ...){
    list(data = object, data_type = 'image')
}

## Assume all point processes in data will be evaluated at the same set of coords
     ## single   



#' @exportS3Method
c.covariate_placeholder <- function(x, y) {
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
    cat('With ',  nrow(attr(object, 'coords')), ' coordinates.')
    cat('\n')
    print(strip_meta(object))
}


#' @export
`[.covariate_placeholder` <- function(object, ...){
  covariate_placeholder(extract_data(object), extract_coords(object))
}

#' @export
extract_data <- function(object) {
  attr(object, 'data')
}

#' @export
extract_coords <- function(object) {
    strip_formatting(object)
}

#' @export
remap <- function(placeholder, new_coords) {
  covariate_placeholder(extract_data(placeholder), new_coords)
}



