# Calculates multiple scores

Calculates multiple scores for every variable pair in a dataset.

## Usage

``` r
pairwise_multi(
  d,
  scores = c("pair_cor", "pair_spearman", "pair_dcor", "pair_mine", "pair_ace",
    "pair_cancor", "pair_nmi", "pair_uncertainty", "pair_chi"),
  handle.na = TRUE,
  warnings = TRUE
)
```

## Arguments

- d:

  dataframe

- scores:

  a character vector naming functions returning a `pairwise` from a
  dataset.

- handle.na:

  If TRUE uses pairwise complete observations to calculate pairwise
  score, otherwise NAs not handled.

- warnings:

  If TRUE, generates a warning for datasets of one row, one column, or
  with constant variables.

## Value

tibble of class "pairwise"

## Examples

``` r
iris1 <- iris
iris1$Sepal.Length <- cut(iris1$Sepal.Length,3)
pairwise_multi(iris1)
#> # A tibble: 44 × 6
#>    x            y            score    group value pair_type
#>    <chr>        <chr>        <chr>    <chr> <dbl> <chr>    
#>  1 Petal.Length Petal.Width  pearson  all   0.963 nn       
#>  2 Petal.Length Petal.Width  spearman all   0.938 nn       
#>  3 Petal.Length Petal.Width  dcor     all   0.974 nn       
#>  4 Petal.Length Petal.Width  MIC      all   0.918 nn       
#>  5 Petal.Length Petal.Width  ace      all   0.989 nn       
#>  6 Petal.Length Petal.Width  cancor   all   0.963 nn       
#>  7 Petal.Length Petal.Width  nmi      all   0.835 nn       
#>  8 Petal.Length Sepal.Length ace      all   0.865 fn       
#>  9 Petal.Length Sepal.Length cancor   all   0.860 fn       
#> 10 Petal.Length Sepal.Length nmi      all   0.539 fn       
#> # ℹ 34 more rows
```
