test_that("Helper file running",
{
    expect_no_error(covar)
    expect_s3_class(covar, 'Pcov')
})


test_that("Smooth construction",
{
    object <- s(covar, bs = 'conv', k = 17)
    pcovar <- covariate_placeholder(covar, test_points)
    data <- list(covar = pcovar)

    knots <- NULL
    smooth <- smooth.construct(object = object, data = data, knots = knots)

    expect_equal(ncol(smooth$X), 17)
    expect_equal( class(smooth$internal_basis), "Bspline2.smooth")
})

test_that("Smooth distance cutoff",
{
    object <- s(covar, bs = 'conv', k = 17, xt = list(max_dist_prop = 0.05 ))
    pcovar <- covariate_placeholder(covar, test_points)
    data <- list(covar = pcovar)

    smooth <- smooth.construct(object = object, data = data, knots = knots)
    
})
          



test_that("Smooth Distance Selection", 
{
    object <- s(covar, bs = 'conv', k = 17, xt = list(max_dist_prop = 0.05 ))
    pcovar <- covariate_placeholder(covar, test_points)
    data <- list(covar = pcovar)
    smooth <- smooth.construct(object = object, data = data, knots = knots)
    
})


test_that("Model setup works",
{
    pcovar <- covariate_placeholder(covar, test_points)

    df <- tibble::tibble(
        y = rnorm(length(pcovar)),
        x = pcovar
    )

        object <- s(covar, bs = 'conv', k = 17, xt = list(max_dist_prop = 0.05 ))
    gam(y ~ s(x, bs = 'conv'), data = df)
    
    
    odata <- list(covar = pcovar)
    smooth <- smooth.construct(object = object, data = data, knots = knots)
    
})
          

          
