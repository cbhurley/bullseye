# Kendall's W for association between ordinal factors.

Calculates Kendall's tau W every factor variable pair in a dataset.

## Usage

``` r
pair_tauW(d, handle.na = TRUE, warnings = TRUE, ...)
```

## Arguments

- d:

  A dataframe

- handle.na:

  ignored. Pairwise complete observations are used automatically.

- warnings:

  If TRUE, generates a warning for datasets of one row, one column, or
  with constant variables.

- ...:

  other arguments

## Value

A tibble of class `pairwise` with factor pairs, or NULL if there are not
at least two factor variables

## Details

Calculated using
[`KendallW`](https://andrisignorell.github.io/DescTools/reference/KendallW.html).
Assumes factor levels are in the given order. NAs are automatically
handled by pairwise omit.

## Examples

``` r
 d <- data.frame(x=rnorm(20), 
                 y=factor(sample(3,20, replace=TRUE)), 
                 z=factor(sample(2,20, replace=TRUE)))
 pair_tauW(d)
#> # A tibble: 1 × 6
#>   x     y     score group value pair_type
#>   <chr> <chr> <chr> <chr> <dbl> <chr>    
#> 1 y     z     tauW  all   0.403 ff       
```
