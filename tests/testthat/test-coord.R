
test_that("norm works",
{
    x <- coord(sqrt(2), sqrt(2))
    expect_equal(Norm(x), 2)
}
)
