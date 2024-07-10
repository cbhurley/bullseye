
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bullseye

<!-- badges: start -->
<!-- badges: end -->

The goal of bullseye is to provide a tidy data structure and
visualisations for multiple or grouped variable correlations, general
association measures and other pairwise scores suitable for numerical,
ordinal and nominal variables.

## Installation

You can install the development version of bullseye from
[GitHub](https://github.com/) with:

``` r
# still in progress
# install.packages("devtools")
# devtools::install_github("cbhurley/bullseye")
```

## Build a `pairwise` data structure

``` r
irisc <- pair_scores(iris, by = "Species") 
irisc
```

This calculates correlations for every level of species in the data.

## Visualise the correlations

``` r
plot_pairwise(irisc)
```
