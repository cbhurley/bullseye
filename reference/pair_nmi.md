# Normalized mutual information

Calculates normalized mutual information for every numeric or factor or
mixed variable pair in a dataset.

## Usage

``` r
pair_nmi(d, handle.na = TRUE, warnings = TRUE, ...)
```

## Arguments

- d:

  A dataframe

- handle.na:

  If TRUE uses pairwise complete observations to calculate normalized
  mutual information, otherwise NAs not handled.

- warnings:

  If TRUE, generates a warning for datasets of one row, one column, or
  with constant variables.

- ...:

  other arguments

## Value

A tibble of class `pairwise`

## Details

The normalized mutual information is calculated using
[`maxNMI`](https://rdrr.io/pkg/linkspotter/man/maxNMI.html) from
linkpotter package

## Examples

``` r
if (requireNamespace("linkspotter", quietly = TRUE)) { 
   pair_nmi(iris)
}
#> # A tibble: 10 × 6
#>    x            y            score group value pair_type
#>    <chr>        <chr>        <chr> <chr> <dbl> <chr>    
#>  1 Petal.Length Sepal.Length nmi   all   0.699 nn       
#>  2 Petal.Width  Sepal.Length nmi   all   0.632 nn       
#>  3 Sepal.Length Sepal.Width  nmi   all   0.203 nn       
#>  4 Petal.Length Sepal.Width  nmi   all   0.379 nn       
#>  5 Petal.Width  Sepal.Width  nmi   all   0.370 nn       
#>  6 Petal.Length Petal.Width  nmi   all   0.835 nn       
#>  7 Sepal.Length Species      nmi   all   0.487 fn       
#>  8 Sepal.Width  Species      nmi   all   0.261 fn       
#>  9 Petal.Length Species      nmi   all   0.870 fn       
#> 10 Petal.Width  Species      nmi   all   0.892 fn       
```
