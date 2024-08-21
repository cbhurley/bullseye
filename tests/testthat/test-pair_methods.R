

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
  if (requireNamespace("energy", quietly = TRUE)){
    p <- pair_dcor(iris[1:4,])
    expect_s3_class(p, "pairwise")
    expect_identical(dim(p), c(6L,6L))
  }
})


test_that("pair mine", {
  if (requireNamespace("minerva", quietly = TRUE)){
    p <- pair_mine(iris[c(1,2,51,52,53),])
    expect_s3_class(p, "pairwise")
    expect_identical(dim(p), c(6L,6L))
  }
})

test_that("pair nmi", {
  if (requireNamespace("linkspotter", quietly = TRUE)){
    p <- pair_nmi(iris[1:4,])
    expect_s3_class(p, "pairwise")
    expect_identical(dim(p), c(10L,6L))
  }
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
  if (requireNamespace("DescTools", quietly = TRUE)){
    iris1 <- iris[c(1,2,53,55),]
    iris1$Sepal.Length <- cut(iris1$Sepal.Length,2)
    p <- pair_tauB(droplevels(iris1))
    expect_s3_class(p, "pairwise")
    expect_identical(dim(p), c(1L,6L))
    p <- pair_tauA(droplevels(iris1))
    expect_s3_class(p, "pairwise")
    expect_identical(dim(p), c(1L,6L))
    p <- pair_tauC(droplevels(iris1))
    expect_s3_class(p, "pairwise")
    expect_identical(dim(p), c(1L,6L))
    p <- pair_tauW(droplevels(iris1))
    expect_s3_class(p, "pairwise")
    expect_identical(dim(p), c(1L,6L))
  }
})


test_that("pair uncertainty", {
  if (requireNamespace("DescTools", quietly = TRUE)){
    iris1 <- iris[c(1,2,53,55),]
    iris1$Sepal.Length <- cut(iris1$Sepal.Length,2)
    p <- pair_uncertainty(droplevels(iris1))
    expect_s3_class(p, "pairwise")
    expect_identical(dim(p), c(1L,6L))
  }
})

test_that("pair gkTau", {
  iris1 <- iris[c(1,2,53,55),]
  iris1$Sepal.Length <- cut(iris1$Sepal.Length,2)
  p <- pair_gkTau(droplevels(iris1))
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(1L,6L))
})


test_that("pair gkGamma", {
  if (requireNamespace("DescTools", quietly = TRUE)){
    iris1 <- iris[c(1,2,53,55),]
    iris1$Sepal.Length <- cut(iris1$Sepal.Length,2)
    p <- pair_gkGamma(droplevels(iris1))
    expect_s3_class(p, "pairwise")
    expect_identical(dim(p), c(1L,6L))
  }
})

test_that("pair chi", {
  if (requireNamespace("DescTools", quietly = TRUE)){
    iris1 <- iris[c(1,2,53,55),]
    iris1$Sepal.Length <- cut(iris1$Sepal.Length,2)
    p <- pair_chi(droplevels(iris1))
    expect_s3_class(p, "pairwise")
    expect_identical(dim(p), c(1L,6L))
  }
})

test_that("pair ace", {
  if (requireNamespace("acepack", quietly = TRUE)){
    p <- pair_ace(iris[c(1,2,53,55),])
    expect_s3_class(p, "pairwise")
    expect_identical(dim(p), c(10L,6L))
  }
})


test_that("pair multi", {
  p <- pairwise_multi(iris)
  expect_s3_class(p, "pairwise")
  expect_identical(dim(p), c(54L,6L))
})

