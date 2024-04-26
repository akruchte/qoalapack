## cc <- tigris::counties(state = 'PA')

## cc <- cc |> select(water = AWATER)


## library(sf)
## library(mgcv)
## library(tidyverse)
## bb <- st_bbox(cc)


gridx <- seq(from = bb$xmin, to = bb$xmax, length.out = 128)
gridy <- seq(from = bb$ymin, to = bb$ymax, length.out = 128)
stepx <- gridx[2] - gridx[1]
stepy <- gridy[2] - gridy[1]

grid <- expand_grid(x = gridx, y = gridy)

points <- st_as_sf(grid, coords = c('x', 'y'))
st_crs(points) <- st_crs(cc)

within <- st_within(points, cc)

included <- lengths(within) > 0

points <- mutate(points[included,], membership = as.integer(within[included]))
coords <- st_coordinates(points)
coords <- as.data.frame(coords)
names(coords) <- c('x', 'y')
internal_basis <- smooth.construct(s(x,y, bs = 'tp'), data = coords, knots = NULL)

X <- internal_basis$X * stepx * stepy

integrated <- do.call(rbind, lapply(split.data.frame(X, points$membership), colSums))
colSums(X[points$membership == 6,])


point_membership <- function(data) {
  stopifnot(names(data) == c('x', 'y', 'membership'))
  structure(data)
}
smooth.construct.int.smooth.spec <- function(object, data, knots){
internal_basis <- smooth.construct(s(x,y, bs = 'tp'), data = coords, knots = NULL)
}

