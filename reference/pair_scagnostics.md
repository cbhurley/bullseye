# Graph-theoretic scagnostics values

Calculates scagnostic values for every numeric variable pair in a
dataset.

## Usage

``` r
pair_scagnostics(
  d,
  scagnostic = c("Outlying", "Skewed", "Clumpy", "Sparse", "Striated", "Convex",
    "Skinny", "Stringy", "Monotonic"),
  handle.na = TRUE,
  warnings = TRUE,
  ...
)
```

## Arguments

- d:

  A dataframe

- scagnostic:

  a character vector for the scagnostic to be calculated. Subset of
  "Outlying", "Stringy", "Striated", "Clumpy", "Sparse", "Skewed",
  "Convex", "Skinny" or "Monotonic"

- handle.na:

  If TRUE uses pairwise complete observations.

- warnings:

  If TRUE, generates a warning for datasets of one row, one column, or
  with constant variables.

- ...:

  other arguments

## Value

A tibble of class `pairwise` with scagnostic values for every numeric
variable pair, or NULL if there are not at least two numeric variables

## Details

The scagnostic values are calculated using
[`scagnostics`](https://rdrr.io/pkg/scagnostics/man/scagnostics.html)
function from the `scagnostics` package.

## References

Wilkinson, Leland, Anushka Anand, and Robert Grossman. "Graph-theoretic
scagnostics." Information Visualization, IEEE Symposium on. IEEE
Computer Society, 2005

## Examples

``` r
pair_scagnostics(iris)
#> # A tibble: 54 × 6
#>    x            y            score     group  value pair_type
#>    <chr>        <chr>        <chr>     <chr>  <dbl> <chr>    
#>  1 Petal.Length Petal.Width  Outlying  all   0.0484 nn       
#>  2 Petal.Length Petal.Width  Skewed    all   0.428  nn       
#>  3 Petal.Length Petal.Width  Clumpy    all   0.549  nn       
#>  4 Petal.Length Petal.Width  Sparse    all   0.0507 nn       
#>  5 Petal.Length Petal.Width  Striated  all   0.108  nn       
#>  6 Petal.Length Petal.Width  Convex    all   0.316  nn       
#>  7 Petal.Length Petal.Width  Skinny    all   0.662  nn       
#>  8 Petal.Length Petal.Width  Stringy   all   0.393  nn       
#>  9 Petal.Length Petal.Width  Monotonic all   0.889  nn       
#> 10 Petal.Length Sepal.Length Outlying  all   0.247  nn       
#> # ℹ 44 more rows
```
