if (FALSE){
ggplot(datp, aes(x, y, fill = res)) + 
    scale_fill_viridis() + 
    geom_tile()

ggplot(datp, aes(x, y, fill = `conv(prod)`)) +
    scale_fill_viridis() + 
    geom_tile()


ggplot(datp, aes(x, y, fill = `conv(for_seroad)`)) +
    scale_fill_viridis() + 
    geom_tile()

ggplot(datp, aes(x, y, fill = `conv(for_prroad)`)) +
    scale_fill_viridis() + 
    geom_tile()





ggsave('pm2.5raster.pdf')

ggplot(datp, aes(x,y, fill = resp)) +
    scale_fill_viridis() + 
    geom_tile()

ggplot(for_plots) +
    scale_fill_viridis() + 
    geom_tile(aes(x,y, fill = residuals)) +
    geom_sf(data = primary_roads)


library(janitor)
countyfiles <- tar_read("countyfiles") |> 
    clean_names() |>
    select(geoid, aland, awater, geometry)


oldxy <- select(dat, x, y)
    
ggplot() +
    geom_sf(data = primary_roads) +
    geom_sf(data = prod_local,  size = 1) +
    geom_sf(data = clairton_coke, color = 'red', size = 2)+
theme_cowplot() +
labs(x = 'Longitude', y = 'Latitude') +
ggtitle("Map of Study Exposure Data")

ggsave('study_data_map.pdf')

ggplot(datp, aes(x,y, fill = pred)) + geom_tile() + scale_fill_viridis()
ggplot(datp, aes(x,y, color = `conv(clairton)`)) + geom_tile()
ggplot(datp, aes(x,y, color = `conv(prod)`)) + geom_tile()

## write_rds(mods, local('air_pollution_mods.rds'))
## write_rds(mods_c, local('air_pollution_mods_c.rds'))
## write_rds(mods_w, local('air_pollution_mods_w.rds'))

prod_local_buffer <- st_buffer(prod_local, dist = as_units(5, 'km'))
prod_local_buffer_10 <- st_buffer(prod_local, dist = as_units(10, 'km'))
prod_local_buffer_25 <- st_buffer(prod_local, dist = as_units(25, 'km'))

comp_points <- lengths(st_within(observation_grid_points, prod_local_buffer)) > 0
comp_points_10 <- lengths(st_within(observation_grid_points, prod_local_buffer_10)) > 0
comp_points_25 <- lengths(st_within(observation_grid_points, prod_local_buffer_25)) > 0
dat_templ$comp_points <- comp_points
dat_templ$comp_points_10 <- comp_points_10
dat_templ$comp_points_25 <- comp_points_25

fp <- for_plots |>  mutate(in_buf = drop_na(dat_templ)$comp_points, in_buf10 = drop_na(dat_templ)$comp_points_10, in_buf25 = drop_na(dat_templ)$comp_points_25)


MM <- model.matrix(model)
library(stringr)
which_cols <- str_detect(colnames(MM), 'prod')
coefs <- coef(model)
vc <- vcov(model)

contrast_vec <- if_else(fp$in_buf, 1/ sum(fp$in_buf), -1 / (sum(!fp$in_buf)))
v <- sqrt(contrast_vec %*% MM %*% vc %*% t(contrast_vec %*% MM %*% vc)) / sqrt(nrow(MM) - length(coefs))
est <- contrast_vec %*% MM %*% coefs


contrast_vec_10 <- if_else(fp$in_buf10, 1/ sum(fp$in_buf10), -1 / (sum(!fp$in_buf10)))
v <- sqrt(contrast_vec_10 %*% MM %*% vc %*% t(contrast_vec_10 %*% MM %*% vc)) / sqrt(nrow(MM) - length(coefs))
est <- contrast_vec_10 %*% MM %*% coefs

contrast_vec_25 <- if_else(fp$in_buf25, 1/ sum(fp$in_buf25), -1 / (sum(!fp$in_buf25)))
v <- sqrt(contrast_vec_25 %*% MM %*% vc %*% t(contrast_vec_25 %*% MM %*% vc)) / sqrt(nrow(MM) - length(coefs))
est <- contrast_vec_25 %*% MM %*% coefs



latlon <- st_transform(prod_local, CRS)

ggplot(for_plots) + geom_raster(aes(x,y, fill = production)) +
theme_cowplot() +
labs(fill = 'Relative PM2.5\n Contribution', x = 'Longitude', 'Latitude') +
scale_fill_viridis() + 
ggtitle('Estimated Effect of Active Production Wells on PM2.5') +
theme(legend.position = 'bottom') +
geom_sf(data = latlon)

ggsave(local('wells_production_plot.pdf'))

}
