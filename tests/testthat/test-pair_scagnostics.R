
test_that("pair scagnostics", {
  skip_if_not(requireNamespace("scagnostics", quietly = TRUE), message = "Package scagnostics not available.")
  p <- pair_scagnostics(iris[c(1,2,53,55),], c("Clumpy", "Striated"))
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(12L,6L))
})