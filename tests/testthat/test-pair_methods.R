

test_that("pair cor", {
  p <- pair_cor(iris[c(1,2,51,52,53),])
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(6L,6L))
})

test_that("pair cancor", {
  p <- pair_cancor(iris[1:4,])
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(10L,6L))
})

test_that("pair dcor", {
  p <- pair_dcor(iris[1:4,])
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(6L,6L))
})


test_that("pair mine", {
  p <- pair_mine(iris[c(1,2,51,52,53),])
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(6L,6L))
})

test_that("pair nmi", {
  p <- pair_nmi(iris[1:4,])
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
  iris1 <- iris[c(1,2,53,55),]
  iris1$Sepal.Length <- cut(iris1$Sepal.Length,2)
  p <- pair_tau(droplevels(iris1))
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(1L,6L))
})


test_that("pair uncertainty", {
  iris1 <- iris[c(1,2,53,55),]
  iris1$Sepal.Length <- cut(iris1$Sepal.Length,2)
  p <- pair_uncertainty(droplevels(iris1))
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(1L,6L))
})

test_that("pair gkTau", {
  iris1 <- iris[c(1,2,53,55),]
  iris1$Sepal.Length <- cut(iris1$Sepal.Length,2)
  p <- pair_gkTau(droplevels(iris1))
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(1L,6L))
})


test_that("pair gkGamma", {
  iris1 <- iris[c(1,2,53,55),]
  iris1$Sepal.Length <- cut(iris1$Sepal.Length,2)
  p <- pair_gkGamma(droplevels(iris1))
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(1L,6L))
})

test_that("pair chi", {
  iris1 <- iris[c(1,2,53,55),]
  iris1$Sepal.Length <- cut(iris1$Sepal.Length,2)
  p <- pair_chi(droplevels(iris1))
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(1L,6L))
})

test_that("pair ace", {
  p <- pair_ace(iris[c(1,2,53,55),])
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(10L,6L))
})


# test_that("pair multi", {
#   p <- pair_multi(iris)
#   expect_s3_class(p, "pairwise")
#   expect_identical(dim(p), c(60L,6L))
# })

