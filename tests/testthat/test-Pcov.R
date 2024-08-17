
test_that("Pcov works",
{
    TD <- spatstat.data::swedishpines
    expect_length(Pcov(TD), 1L)
})

test_that("Pcov works on multiple data",
{
    TD <- spatstat.data::swedishpines
    expect_length(Pcov(TD, TD), 2)
})
