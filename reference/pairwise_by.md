# Constructs a pairwise result for each level of a by variable.

Constructs a pairwise result for each level of a by variable.

## Usage

``` r
pairwise_by(
  d,
  by,
  pair_fun,
  ungrouped = TRUE,
  warnings = TRUE,
  add.nobs = FALSE
)
```

## Arguments

- d:

  a dataframe

- by:

  a character string for the name of the conditioning variable.

- pair_fun:

  A function returning a `pairwise` from a dataset.

- ungrouped:

  If TRUE calculates the ungrouped score in addition to grouped scores.

- warnings:

  If TRUE, generates a warning for datasets of one row, one column, or
  with constant variables.

- add.nobs:

  If TRUE, adds an extra column containing the number of observations
  used for each score.

## Value

tibble of class "pairwise"

## Examples

``` r
pairwise_by(iris, by="Species", pair_cor)
#> # A tibble: 24 × 6
#>    x            y            score   group      value pair_type
#>    <chr>        <chr>        <chr>   <fct>      <dbl> <chr>    
#>  1 Petal.Length Petal.Width  pearson setosa     0.332 nn       
#>  2 Petal.Length Petal.Width  pearson versicolor 0.787 nn       
#>  3 Petal.Length Petal.Width  pearson virginica  0.322 nn       
#>  4 Petal.Length Petal.Width  pearson all        0.963 nn       
#>  5 Petal.Length Sepal.Length pearson setosa     0.267 nn       
#>  6 Petal.Length Sepal.Length pearson versicolor 0.754 nn       
#>  7 Petal.Length Sepal.Length pearson virginica  0.864 nn       
#>  8 Petal.Length Sepal.Length pearson all        0.872 nn       
#>  9 Petal.Length Sepal.Width  pearson setosa     0.178 nn       
#> 10 Petal.Length Sepal.Width  pearson versicolor 0.561 nn       
#> # ℹ 14 more rows

```
