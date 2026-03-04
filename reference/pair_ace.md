# Alternating conditional expectations correlation

Calculates the maximal correlation coefficient from alternating
conditional expectations algorithm for every variable pair in a dataset.

## Usage

``` r
pair_ace(d, handle.na = TRUE, warnings = TRUE, ...)
```

## Arguments

- d:

  A dataframe

- handle.na:

  If TRUE uses pairwise complete observations, otherwise NAs not
  handled.

- warnings:

  If TRUE, generates a warning for datasets of one row, one column, or
  with constant variables.

- ...:

  other arguments

## Value

A tibble of class `pairwise` with a maximal correlation from the
alternating conditional expectations algorithm for every variable pair

## Details

The maximal correlation is calculated using alternating conditional
expectations algorithm which find the transformations of variables such
that the squared correlation is maximised. The
[`ace`](https://rdrr.io/pkg/acepack/man/ace.html) function from
`acepack` package is used for the calculation.

## References

Breiman, Leo, and Jerome H. Friedman. "Estimating optimal
transformations for multiple regression and correlation." Journal of the
American statistical Association 80.391 (1985): 580-598.

## Examples

``` r
pair_ace(iris)
#> # A tibble: 10 × 6
#>    x            y            score group value pair_type
#>    <chr>        <chr>        <chr> <chr> <dbl> <chr>    
#>  1 Petal.Length Sepal.Length ace   all   0.913 nn       
#>  2 Petal.Width  Sepal.Length ace   all   0.865 nn       
#>  3 Sepal.Length Sepal.Width  ace   all   0.584 nn       
#>  4 Petal.Length Sepal.Width  ace   all   0.706 nn       
#>  5 Petal.Width  Sepal.Width  ace   all   0.731 nn       
#>  6 Petal.Length Petal.Width  ace   all   0.989 nn       
#>  7 Sepal.Length Species      ace   all   0.838 fn       
#>  8 Sepal.Width  Species      ace   all   0.679 fn       
#>  9 Petal.Length Species      ace   all   0.994 fn       
#> 10 Petal.Width  Species      ace   all   0.994 fn       
```
