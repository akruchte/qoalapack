
test_that("Pcov works",
{
    TD <- covar
    expect_length(TD, 1L)
})

test_that("Pcov works on multiple data",
{
    TD <- covar
    expect_length(c(TD, TD), 2)
})
