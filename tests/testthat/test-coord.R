x <- coord(sqrt(2), sqrt(2))
test_that("norm works",
          {
              expect_equal(Norm(x), 2)
          }
)
