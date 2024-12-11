#' Data extent
#'
#' Slightly more general than the observation window, this bounds possible regions
#' that a covariate may be evaluable in.
#' @export
extent <- function (ob, ...) {
    UseMethod('extent')
}


