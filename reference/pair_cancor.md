# Canonical correlation

Calculates canonical correlation for every variable pair in a dataset.

## Usage

``` r
pair_cancor(d, handle.na = TRUE, warnings = TRUE, ...)
```

## Arguments

- d:

  A dataframe

- handle.na:

  If TRUE uses pairwise complete observations to calculate correlation
  coefficient,, otherwise NAs not handled.

- warnings:

  If TRUE, generates a warning for datasets of one row, one column, or
  with constant variables.

- ...:

  other arguments

## Value

A tibble of class `pairwise` with canonical correlation for every
numeric or factor or mixed variable pair

## Examples

``` r
pair_cancor(iris)
#> # A tibble: 10 × 6
#>    x            y            score  group value pair_type
#>    <chr>        <chr>        <chr>  <chr> <dbl> <chr>    
#>  1 Petal.Length Sepal.Length cancor all   0.872 nn       
#>  2 Petal.Width  Sepal.Length cancor all   0.818 nn       
#>  3 Sepal.Length Sepal.Width  cancor all   0.118 nn       
#>  4 Petal.Length Sepal.Width  cancor all   0.428 nn       
#>  5 Petal.Width  Sepal.Width  cancor all   0.366 nn       
#>  6 Petal.Length Petal.Width  cancor all   0.963 nn       
#>  7 Sepal.Length Species      cancor all   0.787 fn       
#>  8 Sepal.Width  Species      cancor all   0.633 fn       
#>  9 Petal.Length Species      cancor all   0.970 fn       
#> 10 Petal.Width  Species      cancor all   0.964 fn       
```
