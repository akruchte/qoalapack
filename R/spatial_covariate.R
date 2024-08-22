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


#' @exportS3Method
format.spatial_covariate <- function(object, ...){
    type <- class(object)[[1]]
    rep(type, vctrs::vec_size(object))
}


strip_meta <- c


default_spatial_points_config <- list(npoints = 100L, F2 = 10)


#' Generate Default Evaluation Points
#'
#' This function generates a set of evaluation points at a default set of locations.
#' @param An object of type extent or with an extent method
default_points <-function(ext) {
    SobolSequence::sobolSequence.points(
                       2,
                       default_spatial_points_config$F2,
                       default_spatial_points_config$npoints
                   )
}
