library(fmesher)
library(sf)



points <- cbind(runif(100), runif(100))
dom <- cbind(c(0, 1, 1, 0), c(0, 0, 1, 1))
mesh <- fm_mesh_2d(points, dom)
