


test_that("pairwise scores", {
  p <- pairwise_scores(iris[c(1,2,53,55),])
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(10L,6L))
  p <- pairwise_scores(iris, by="Species")
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(24L,6L))
})