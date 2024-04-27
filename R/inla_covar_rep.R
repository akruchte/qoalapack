

  library(readr)
  library(sf)
  map <- read_rds('R/simulation_utilities/allegheny_region_tracts.rds')

## sample code for generating the fem object
library(tidyverse)
topo <- map |> st_make_valid() |> 
  group_split(id) |> map(st_union) |> 
  map_dfr(st_as_sf)

t1 <- st_as_sfc(topo[1,])

## pa coords--40.4406° N, 79.9959° W
t1 <- st_transform(t1, "EPSG:3857") |> st_simplify(dTolerance = 200)

library(fmesher)
# create mesh
edge <- st_buffer(t1, 3000)
me <- fm_mesh_2d(boundary = edge, interior = t1, max.edge = 2500)
library(INLA)

Amat <- inla.spde.make.A(me)
#fm_basis

spde_mod <- inla.spde2.matern(me, alpha = 2)





#compositional covars


                    