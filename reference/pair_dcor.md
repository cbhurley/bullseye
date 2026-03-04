# Distance correlation

Calculates distance correlation for every numeric variable pair in a
dataset.

## Usage

``` r
pair_dcor(d, handle.na = TRUE, warnings = TRUE, ...)
```

## Arguments

- d:

  A dataframe

- handle.na:

  If TRUE uses pairwise complete observations to calculate distance
  correlation, otherwise NAs not handled.

- warnings:

  If TRUE, generates a warning for datasets of one row, one column, or
  with constant variables.

- ...:

  other arguments

## Value

A tibble of class `pairwise` with distance correlation for every numeric
variable pair, or NULL if there are not at least two numeric variables

## Details

The distance correlation is calculated using
[`dcor2d`](https://rdrr.io/pkg/energy/man/dcov2d.html) from `energy`
package

## Examples

``` r
pair_dcor(iris)
#> # A tibble: 6 × 6
#>   x            y            score group value pair_type
#>   <chr>        <chr>        <chr> <chr> <dbl> <chr>    
#> 1 Petal.Length Sepal.Length dcor  all   0.859 nn       
#> 2 Petal.Width  Sepal.Length dcor  all   0.827 nn       
#> 3 Sepal.Length Sepal.Width  dcor  all   0.311 nn       
#> 4 Petal.Length Sepal.Width  dcor  all   0.542 nn       
#> 5 Petal.Width  Sepal.Width  dcor  all   0.513 nn       
#> 6 Petal.Length Petal.Width  dcor  all   0.974 nn       
```
