
## TODO
default_spatial_points_config <- list(npoints = 100L, F2 = 10)


default_points <-function() {
    SobolSequence::sobolSequence.points(
                       2,
                       default_spatial_points_config$F2,
                       default_spatial_points_config$npoints
                   )
}
