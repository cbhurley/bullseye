# Kendall's correlation

Calculates Kendall's correlation for every numeric variable pair in a
dataset.

## Usage

``` r
pair_kendall(d, handle.na = TRUE, warnings = TRUE, ...)
```

## Arguments

- d:

  A dataframe

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
pair_kendall(iris)
#> # A tibble: 6 × 6
#>   x            y            score   group   value pair_type
#>   <chr>        <chr>        <chr>   <chr>   <dbl> <chr>    
#> 1 Petal.Length Sepal.Length kendall all    0.719  nn       
#> 2 Petal.Width  Sepal.Length kendall all    0.655  nn       
#> 3 Sepal.Length Sepal.Width  kendall all   -0.0770 nn       
#> 4 Petal.Length Sepal.Width  kendall all   -0.186  nn       
#> 5 Petal.Width  Sepal.Width  kendall all   -0.157  nn       
#> 6 Petal.Length Petal.Width  kendall all    0.807  nn       
# same as
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
```
