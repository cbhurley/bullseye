# Pearson, Spearman or Kendall correlation

Calculates one of either pearson, spearman or kendall correlation for
every numeric variable pair in a dataset.

## Usage

``` r
pair_cor(d, method = "pearson", handle.na = TRUE, warnings = TRUE, ...)
```

## Arguments

- d:

  A dataframe

- method:

  A character string for the correlation coefficient to be calculated.
  Either "pearson" (default), "spearman", or "kendall". If the value is
  "all", then all three correlations are calculated.

- handle.na:

  If TRUE uses pairwise complete observations to calculate correlation
  coefficient, otherwise NAs not handled.

- warnings:

  If TRUE, generates a warning for datasets of one row, one column, or
  with constant variables.

- ...:

  other arguments

## Value

A tibble of class `pairwise` with calculated association value for every
numeric variable pair, or NULL if there are not at least two numeric
variables

## See also

See
[`pair_methods`](https://cbhurley.github.io/bullseye/reference/pair_methods.md)
for other score options.

## Examples

``` r
pair_cor(iris)
#> # A tibble: 6 × 6
#>   x            y            score   group  value pair_type
#>   <chr>        <chr>        <chr>   <chr>  <dbl> <chr>    
#> 1 Petal.Length Sepal.Length pearson all    0.872 nn       
#> 2 Petal.Width  Sepal.Length pearson all    0.818 nn       
#> 3 Sepal.Length Sepal.Width  pearson all   -0.118 nn       
#> 4 Petal.Length Sepal.Width  pearson all   -0.428 nn       
#> 5 Petal.Width  Sepal.Width  pearson all   -0.366 nn       
#> 6 Petal.Length Petal.Width  pearson all    0.963 nn       
pair_cor(iris, method="kendall")
#> # A tibble: 6 × 6
#>   x            y            score   group   value pair_type
#>   <chr>        <chr>        <chr>   <chr>   <dbl> <chr>    
#> 1 Petal.Length Sepal.Length kendall all    0.719  nn       
#> 2 Petal.Width  Sepal.Length kendall all    0.655  nn       
#> 3 Sepal.Length Sepal.Width  kendall all   -0.0770 nn       
#> 4 Petal.Length Sepal.Width  kendall all   -0.186  nn       
#> 5 Petal.Width  Sepal.Width  kendall all   -0.157  nn       
#> 6 Petal.Length Petal.Width  kendall all    0.807  nn       
pair_cor(iris, method="spearman")
#> # A tibble: 6 × 6
#>   x            y            score    group  value pair_type
#>   <chr>        <chr>        <chr>    <chr>  <dbl> <chr>    
#> 1 Petal.Length Sepal.Length spearman all    0.882 nn       
#> 2 Petal.Width  Sepal.Length spearman all    0.834 nn       
#> 3 Sepal.Length Sepal.Width  spearman all   -0.167 nn       
#> 4 Petal.Length Sepal.Width  spearman all   -0.310 nn       
#> 5 Petal.Width  Sepal.Width  spearman all   -0.289 nn       
#> 6 Petal.Length Petal.Width  spearman all    0.938 nn       
pair_cor(iris, method="all")
#> # A tibble: 18 × 6
#>    x            y            score    group   value pair_type
#>    <chr>        <chr>        <chr>    <chr>   <dbl> <chr>    
#>  1 Petal.Length Sepal.Length pearson  all    0.872  nn       
#>  2 Petal.Width  Sepal.Length pearson  all    0.818  nn       
#>  3 Sepal.Length Sepal.Width  pearson  all   -0.118  nn       
#>  4 Petal.Length Sepal.Width  pearson  all   -0.428  nn       
#>  5 Petal.Width  Sepal.Width  pearson  all   -0.366  nn       
#>  6 Petal.Length Petal.Width  pearson  all    0.963  nn       
#>  7 Petal.Length Sepal.Length spearman all    0.882  nn       
#>  8 Petal.Width  Sepal.Length spearman all    0.834  nn       
#>  9 Sepal.Length Sepal.Width  spearman all   -0.167  nn       
#> 10 Petal.Length Sepal.Width  spearman all   -0.310  nn       
#> 11 Petal.Width  Sepal.Width  spearman all   -0.289  nn       
#> 12 Petal.Length Petal.Width  spearman all    0.938  nn       
#> 13 Petal.Length Sepal.Length kendall  all    0.719  nn       
#> 14 Petal.Width  Sepal.Length kendall  all    0.655  nn       
#> 15 Sepal.Length Sepal.Width  kendall  all   -0.0770 nn       
#> 16 Petal.Length Sepal.Width  kendall  all   -0.186  nn       
#> 17 Petal.Width  Sepal.Width  kendall  all   -0.157  nn       
#> 18 Petal.Length Petal.Width  kendall  all    0.807  nn       
```
