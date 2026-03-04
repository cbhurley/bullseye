# Polyserial correlation

Calculates Polyserial correlation using for every factor-numeric
variable pair in a dataset.

## Usage

``` r
pair_polyserial(d, handle.na = TRUE, warnings = TRUE, ...)
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

A tibble of class `pairwise` with polyserial correlation for
factor-numeric pairs, or NULL if there are not at least one such pair.

## Details

The polyserial correlation is calculated using the
[`polyserial`](https://rdrr.io/pkg/polycor/man/polyserial.html) function
from the `polycor` package, and assumes factor levels are in the given
order. NAs are automatically handled by pairwise omit.

## Examples

``` r
pair_polyserial(iris)
#> # A tibble: 4 × 6
#>   x            y       score      group  value pair_type
#>   <chr>        <chr>   <chr>      <chr>  <dbl> <chr>    
#> 1 Sepal.Length Species polyserial all    0.879 fn       
#> 2 Sepal.Width  Species polyserial all   -0.479 fn       
#> 3 Petal.Length Species polyserial all    1.000 fn       
#> 4 Petal.Width  Species polyserial all    1.000 fn       
```
