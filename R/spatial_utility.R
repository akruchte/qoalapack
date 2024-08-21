
#' @export
distance_raster <- function(dims, xmax, xmin, ymin, ymax){

    xc <- (xmax - xmin) / 2.0
    yc <- (ymax - ymin) / 2.0

    xgrid <- seq(from = xmin - xc, to = xmax - xc, length.out = dims[1])
    ygrid <- seq(from = ymin - yc, to = ymax - yc, length.out = dims[2])
    dists <- outer(xgrid, ygrid,
                   \(x,y) sqrt(x^2 + y^2))
}

