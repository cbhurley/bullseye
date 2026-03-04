# A generic function to create a data structure for summarising variable pairs in a dataset

Creates a data structure for every variable pair in a dataset.

## Usage

``` r
pairwise(x, score = NA_character_, pair_type = NA_character_)

# S3 method for class 'matrix'
pairwise(x, score = NA_character_, pair_type = NA_character_)

# S3 method for class 'data.frame'
pairwise(x, score = NA_character_, pair_type = NA_character_)

# S3 method for class 'easycorrelation'
pairwise(x, score = NA_character_, pair_type = NA_character_)

as.pairwise(x, score = NA_character_, pair_type = NA_character_)
```

## Arguments

- x:

  A dataframe or symmetric matrix.

- score:

  a character string indicating the value of association

- pair_type:

  a character string specifying the type of variable pair, should be
  either "nn", "fn", "ff", for a numeric-numeric pair, factor-numeric
  pair, or factor-factor pair, or NA if unknown.

## Value

A tbl_df of class `pairwise` for pairs of variables with a column
`value` for the score value, `score` for a type of association value and
`pair_type` for the type of variable pair.

## Details

The `pairwise` class has columns x and y for (ordered pairs) of
variables, where x \< y. The column score has the name of the summary
measure used for the two variables, and the column value has the
associated value. The group column defaults to "all", meaning summary
measures apply to the complete dataset, otherwise it describes a subset
of the data. The functions `pair_*` calculate pairwise tibbles for the
summary measure named by `*`, eg
[`pair_cor()`](https://cbhurley.github.io/bullseye/reference/pair_cor.md),
[`pair_cancor()`](https://cbhurley.github.io/bullseye/reference/pair_cancor.md).
The functions
[`pairwise_scores()`](https://cbhurley.github.io/bullseye/reference/pairwise_scores.md)
and
[`pairwise_by()`](https://cbhurley.github.io/bullseye/reference/pairwise_by.md)
calculate pairwise tibbles for levels of a grouping variable. The
function
[`pairwise_multi()`](https://cbhurley.github.io/bullseye/reference/pairwise_multi.md)
calculates a pairwise_tibble for multiple named scores. The pairwise
tibble has at most one row for each combination of x, y, score and
group. This is checked prior to plotting by `plot.pairwise`. Note that
the pair_type column is included for information purposes, but it is not
currently used by `plot.pairwise`.

## Methods (by class)

- `pairwise(matrix)`: pairwise method

- `pairwise(data.frame)`: pairwise method

- `pairwise(easycorrelation)`: pairwise method

## Functions

- `as.pairwise()`: Same as `pairwise`

## Examples

``` r
pairwise(cor(iris[,1:4]), score="pearson")
#> # A tibble: 6 × 6
#>   x            y            score   group  value pair_type
#>   <chr>        <chr>        <chr>   <chr>  <dbl> <chr>    
#> 1 Petal.Length Sepal.Length pearson all    0.872 NA       
#> 2 Petal.Width  Sepal.Length pearson all    0.818 NA       
#> 3 Sepal.Length Sepal.Width  pearson all   -0.118 NA       
#> 4 Petal.Length Sepal.Width  pearson all   -0.428 NA       
#> 5 Petal.Width  Sepal.Width  pearson all   -0.366 NA       
#> 6 Petal.Length Petal.Width  pearson all    0.963 NA       
pairwise(iris)
#> # A tibble: 10 × 6
#>    x            y            score group value pair_type
#>    <chr>        <chr>        <chr> <chr> <dbl> <chr>    
#>  1 Petal.Length Sepal.Length NA    all      NA nn       
#>  2 Petal.Width  Sepal.Length NA    all      NA nn       
#>  3 Sepal.Length Sepal.Width  NA    all      NA nn       
#>  4 Petal.Length Sepal.Width  NA    all      NA nn       
#>  5 Petal.Width  Sepal.Width  NA    all      NA nn       
#>  6 Petal.Length Petal.Width  NA    all      NA nn       
#>  7 Sepal.Length Species      NA    all      NA fn       
#>  8 Sepal.Width  Species      NA    all      NA fn       
#>  9 Petal.Length Species      NA    all      NA fn       
#> 10 Petal.Width  Species      NA    all      NA fn       
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
pairwise_scores(iris, by="Species")
#> # A tibble: 24 × 6
#>    x            y            score   group      value pair_type
#>    <chr>        <chr>        <chr>   <fct>      <dbl> <chr>    
#>  1 Petal.Length Sepal.Length pearson setosa     0.267 nn       
#>  2 Petal.Width  Sepal.Length pearson setosa     0.278 nn       
#>  3 Sepal.Length Sepal.Width  pearson setosa     0.743 nn       
#>  4 Petal.Length Sepal.Width  pearson setosa     0.178 nn       
#>  5 Petal.Width  Sepal.Width  pearson setosa     0.233 nn       
#>  6 Petal.Length Petal.Width  pearson setosa     0.332 nn       
#>  7 Petal.Length Sepal.Length pearson versicolor 0.754 nn       
#>  8 Petal.Width  Sepal.Length pearson versicolor 0.546 nn       
#>  9 Sepal.Length Sepal.Width  pearson versicolor 0.526 nn       
#> 10 Petal.Length Sepal.Width  pearson versicolor 0.561 nn       
#> # ℹ 14 more rows
pairwise_multi(iris)
#> # A tibble: 54 × 6
#>    x            y            score    group value pair_type
#>    <chr>        <chr>        <chr>    <chr> <dbl> <chr>    
#>  1 Petal.Length Petal.Width  pearson  all   0.963 nn       
#>  2 Petal.Length Petal.Width  spearman all   0.938 nn       
#>  3 Petal.Length Petal.Width  dcor     all   0.974 nn       
#>  4 Petal.Length Petal.Width  MIC      all   0.918 nn       
#>  5 Petal.Length Petal.Width  ace      all   0.989 nn       
#>  6 Petal.Length Petal.Width  cancor   all   0.963 nn       
#>  7 Petal.Length Petal.Width  nmi      all   0.835 nn       
#>  8 Petal.Length Sepal.Length pearson  all   0.872 nn       
#>  9 Petal.Length Sepal.Length spearman all   0.882 nn       
#> 10 Petal.Length Sepal.Length dcor     all   0.859 nn       
#> # ℹ 44 more rows
```
