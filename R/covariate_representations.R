

## area to point representation

## a2akriging as a canned option.

## centroid based kernel weighting

## gam based model smooth continous surface

## mock points

library(spatstat)
library(sf)
library(compositions)

point_data <- st_as_sf(swedishpines)[-1,]

vor <- st_voronoi(st_combine(point_data))
vorpoly <- st_collection_extract(vor, 'POLYGON')

centroid_coords <- st_coordinates(st_centroid(vorpoly))

demo1 <- function(x,y) 2 * (x - 0.5) ^2 - 0.3 * x + 0.1 * y
demo2 <- function(x,y) 1.1 * (x - 0.5) ^2 - 0.3 * x + 0.07 * y

demos <- mutate(as_tibble(centroid_coords), demo1 = demo1(X, Y), demo2 = demo2(X,Y)) |>
    mutate(demos = ilrInv(cbind(demo1, demo2))) |>
    mutate(demos = as.data.frame(demos)) |>
    unpack(demos) |>
    janitor::clean_names()

vorpoly <- st_sf(region = vorpoly)
regions <- vorpoly |> bind_cols(demos)


## x1 as default category
## perform multinom fractional logit seperately by categories


## gam(x31 ~ s(x,y), data = mutate(demos, x31 = x3 / (x3 + x1)), family = binomial())

## gam(x21 ~ s(x,y), data = mutate(demos, x21 = x3 / (x2 + x1)), family = binomial())



    




## vorpoly



                       


## various types of deconvolution based methods
