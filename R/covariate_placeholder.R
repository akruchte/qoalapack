## covariate placeholders should carry the relevant information regarding
## the appropriate single level entity  information needed in model fitting (such as mgcv::gam)
## In this case that will most likely be a quadrature point representation
## cases include an age, sex, geo-coordinate (x,y), or possibly higher order coordinates, (x,y,t, w) for extra w
## These should all be based on vctrs and preserve the relevant attributes needed for model setup


## currently hacky implementation that should be improved
## also doesn't seem to be working anyways
## the goal of these is to get mgcv to work in a friendly way
## min.covariate_placeholder <- function(pl, ...) {min(extract_coords(pl)$x, ...)}
## max.covariate_placeholder <- function(pl, ...) {max(extract_coords(pl)$x, ...)}

placeholder_value <- function(value) {
    structure(value, class = c('placeholder_value', 'numeric'))
}
print.placeholder_value <- function(value) {
    cat(str_glue('(({value}))\n\n'))
}

covariate_placeholder <- function(data, coords) {
    structure(rep(placeholder_value(1), nrow(coords)),
              class = c('covariate_placeholder', 'numeric'),
              coords = coords,
              data = data)
}


print.covariate_placeholder <- function(object){
    cat('A Covariate Placeholder\n')
}
`[.covariate_placeholder` <- function(object, ...){
    covariate_placeholder(extract_data(object), extract_coords(object)[...,])
}


extract_data <- function(object) {
    attr(object, 'data')
}

extract_coords <- function(object) {
    attr(object, 'coords')
}

remap <- function(placeholder, new_coords) {
    covariate_placeholder(extract_data(placeholder), new_coords)
}

