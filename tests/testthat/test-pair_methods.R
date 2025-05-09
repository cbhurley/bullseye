
test_that("pair cor", {
  p <- pair_cor(iris[c(1,2,51,52,53),])
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(6L,6L))
})

test_that("pair cancor", {
  p <- pair_cancor(iris[c(1,2,51,52,53),])
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(10L,6L))
})

test_that("pair dcor", {
  skip_if_not(requireNamespace("energy", quietly = TRUE), message = "Package energy not available.")
    p <- pair_dcor(iris[c(1,2,51,52,53),])
    expect_s3_class(p, "pairwise")
    expect_identical(dim(p), c(6L,6L))
})


test_that("pair mine", {
  skip_if_not(requireNamespace("minerva", quietly = TRUE), message = "Package minerva not available.")
  p <- pair_mine(iris[c(1,2,51,52,53),])
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(6L,6L))
})

test_that("pair nmi", {
  skip_if_not(requireNamespace("linkspotter", quietly = TRUE), message = "Package linkspotter not available.")
  p <- pair_nmi(iris[c(1,2,51,52,53),])
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(10L,6L))
})

test_that("pair polychor", {
  iris1 <- iris[c(1,2,53,55),]
  iris1$Sepal.Length <- cut(iris1$Sepal.Length,2)
  p <- pair_polychor(droplevels(iris1))
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(1L,6L))
})

test_that("pair polyserial", {
  iris1 <- iris[c(1,2,53,55),]
  iris1$Sepal.Length <- cut(iris1$Sepal.Length,2)
  p <- pair_polyserial(droplevels(iris1))
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(6L,6L))
})

test_that("pair tau", {
  skip_if_not(requireNamespace("DescTools", quietly = TRUE), message = "Package DescTools not available.")
  iris1 <- iris[c(1,2,53,55),]
  iris1$Sepal.Length <- cut(iris1$Sepal.Length,2)
  iris1 <- droplevels(iris1)
  p <- pair_tauB(iris1)
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(1L,6L))
  p <- pair_tauA(iris1)
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(1L,6L))
  p <- pair_tauC(iris1)
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(1L,6L))
  p <- pair_tauW(iris1)
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(1L,6L))
})


test_that("pair uncertainty", {
  skip_if_not(requireNamespace("DescTools", quietly = TRUE), message = "Package DescTools not available.")
  iris1 <- iris[c(1,2,53,55),]
  iris1$Sepal.Length <- cut(iris1$Sepal.Length,2)
  p <- pair_uncertainty(droplevels(iris1))
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(1L,6L))
})

test_that("pair gkTau", {
  skip_if_not(requireNamespace("DescTools", quietly = TRUE), message = "Package DescTools not available.")
  iris1 <- iris[c(1,2,53,55),]
  iris1$Sepal.Length <- cut(iris1$Sepal.Length,2)
  p <- pair_gkTau(droplevels(iris1))
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(1L,6L))
})


test_that("pair gkGamma", {
  skip_if_not(requireNamespace("DescTools", quietly = TRUE), message = "Package DescTools not available.")
  iris1 <- iris[c(1,2,53,55),]
  iris1$Sepal.Length <- cut(iris1$Sepal.Length,2)
  p <- pair_gkGamma(droplevels(iris1))
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(1L,6L))
})

test_that("pair chi", {
  skip_if_not(requireNamespace("DescTools", quietly = TRUE), message = "Package DescTools not available.")
  iris1 <- iris[c(1,2,53,55),]
  iris1$Sepal.Length <- cut(iris1$Sepal.Length,2)
  p <- pair_chi(droplevels(iris1))
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(1L,6L))
})

test_that("pair ace", {
  skip_if_not(requireNamespace("acepack", quietly = TRUE), message = "Package acepack not available.")
  p <- pair_ace(iris)
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(10L,6L))
})


test_that("pair multi", {
  p <- pairwise_multi(iris, scores=c("pair_cor", "pair_cancor"))
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(16L,6L))
})

