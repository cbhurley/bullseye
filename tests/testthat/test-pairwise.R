

test_that("pairwise", {
  p <- pairwise(iris[,1:4])
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(6L,6L))
})

test_that("pairwise from matrix", {
  p <- pairwise(cor(iris[,1:4]))
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(6L,6L))
})
