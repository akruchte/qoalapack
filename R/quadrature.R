library(sf)
## construct quadrature schemes for point process valued outcomes
## quadscheme from spatstat only handles the purely 2d case, we need to handle temporal, but also other mark variables such as age, sex, etc


## crude first method for now
## directly use spatstat on the unioned points
quadrature <- function(outcomes, ...){

    
    quad_points <- outcomes
    ## quad_grid <- 
    
}

## Ideally this should be able to be constructed adaptively based on data sparsity etc


## window should be provided as a POLYGON or MULTIPOLYGON object
## nx,ny are the (approximate) number of points in the x and y dimensions respectively
## returns a grid of sf POINT objects representing control points for quadrature scheme
control_grid <- function(window, nx, ny, ...){
    st_make_grid(window, n = c(nx, ny), what = 'centers')
}

control_points <- function(window, n, ...){
    st_sample(window, n)
}



