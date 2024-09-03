library(mgcv)

covar <- Pcov(spatstat.data::swedishpines)

window <- spatstat.data::swedishpines$window
set.seed(42)
eval_points <- spatstat.random::runifpoint(100, win = window)
test_points <- coord(spatstat.geom::coords(eval_points))
