## code to prepare `DATASET` dataset goes here
library(tidyverse)
library(tigris)
library(sf)

allegheny <- counties(state = 'PA') |> filter(NAME == 'Allegheny')
tracts <- tracts(state = 'PA', county = 'Allegheny')
majorwater <- area_water(state = 'PA', county = 'Allegheny')


water <- st_union(majorwater)

allegheny <- st_union(allegheny)
allegheny <- st_difference(allegheny, water)

allpolies <- st_cast(allegheny, 'POLYGON')
major_regions <- st_sf(geometry = allpolies) |>
    mutate(area = st_area(geometry)) |>
    arrange(desc(area)) |>
    slice(1:4)

remove_holes <- function(geom) st_polygon(geom[1])
major_regions <- mutate(major_regions, geometry = modify(geometry, remove_holes))
major_regions <- major_regions |> mutate(id = 1:4)
tracts <- select(tracts)
region_tracts <- st_intersection(tracts, major_regions)


usethis::use_data(select(region_tracts, id), overwrite = TRUE)
