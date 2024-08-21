
test_that("pair scagnostics", {
  if (requireNamespace("scagnostics", quietly = TRUE)){
    p <- pair_scagnostics(iris[c(1,2,53,55),], c("Clumpy", "Striated"))
    expect_s3_class(p, "pairwise")
    expect_identical(dim(p), c(12L,6L))
  }
})