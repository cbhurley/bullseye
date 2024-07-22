

test_that("plot pairwise", {
  p <- plot_pairwise(pairwise(iris[1:4,]))
  expect_s3_class(p, "ggplot")
})

test_that("plot pairwise linear", {
  p <- plot_pairwise_linear(pairwise(iris[1:4,]))
  expect_s3_class(p, "ggplot")
})