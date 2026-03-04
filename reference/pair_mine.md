# MINE family values

Calculates MINE family values for every numeric variable pair in a
dataset.

## Usage

``` r
pair_mine(d, method = "MIC", handle.na = TRUE, warnings = TRUE, ...)
```

## Arguments

- d:

  A dataframe

- method:

  character vector for the MINE value to be calculated. Subset of
  "MIC","MAS","MEV","MCN","MICR2", "GMIC", "TIC"

- handle.na:

  If TRUE uses pairwise complete observations to calculate score,
  otherwise NAs not handled.

- warnings:

  If TRUE, generates a warning for datasets of one row, one column, or
  with constant variables.

- ...:

  other arguments

## Value

A tibble of class `pairwise` with scores for numeric variable pairs, or
NULL if there are not at least two numeric variables

## Details

The values are calculated using
[`mine`](https://rdrr.io/pkg/minerva/man/mine.html) from `minerva`

## References

Reshef, David N., et al. "Detecting novel associations in large data
sets." science 334.6062 (2011): 1518-1524

## Examples

``` r
 pair_mine(iris)
#> # A tibble: 6 × 6
#>   x            y            score group value pair_type
#>   <chr>        <chr>        <chr> <chr> <dbl> <chr>    
#> 1 Petal.Length Petal.Width  MIC   all   0.918 nn       
#> 2 Petal.Length Sepal.Length MIC   all   0.768 nn       
#> 3 Petal.Length Sepal.Width  MIC   all   0.439 nn       
#> 4 Petal.Width  Sepal.Length MIC   all   0.668 nn       
#> 5 Petal.Width  Sepal.Width  MIC   all   0.435 nn       
#> 6 Sepal.Length Sepal.Width  MIC   all   0.277 nn       
 pair_mine(iris, method="MAS")
#> # A tibble: 6 × 6
#>   x            y            score group  value pair_type
#>   <chr>        <chr>        <chr> <chr>  <dbl> <chr>    
#> 1 Petal.Length Petal.Width  MAS   all   0.0803 nn       
#> 2 Petal.Length Sepal.Length MAS   all   0.0665 nn       
#> 3 Petal.Length Sepal.Width  MAS   all   0.109  nn       
#> 4 Petal.Width  Sepal.Length MAS   all   0.0365 nn       
#> 5 Petal.Width  Sepal.Width  MAS   all   0.0729 nn       
#> 6 Sepal.Length Sepal.Width  MAS   all   0.0501 nn       
```
