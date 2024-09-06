
#' Direct Convolution
#' @param exposure_points a matrix of exposure point coordinates
#' @param evaluation_points a matrix of evaluation point coordinates
#' @return a m x n array of distances where n is the number of exposure points and m the number of evaluation points
#' @export
direct_conv <- function(exposure_points, evaluation_points) {
    n <- nrow(exposure_points)
    m <- nrow(evaluation_points)

    dists <- array(0, dim = c(m, n))

    for (i in 1:m) {
        for (j in 1:n) {
            dists[i, j] <- sqrt(sum((exposure_points[j,] -  evaluation_points[i,])^2))
        }
    }
    dists
}

