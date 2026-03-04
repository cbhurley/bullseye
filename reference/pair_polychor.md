# Polychoric correlation

Calculates Polychoric correlation using for every factor variable pair
in a dataset.

## Usage

``` r
pair_polychor(d, handle.na = TRUE, warnings = TRUE, ...)
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

A tibble of class `pairwise` with polychoric correlation for factor
pairs, or NULL if there are not at least two factor variables

## Details

The polychoric correlation is calculated using the
[`polychor`](https://rdrr.io/pkg/polycor/man/polychor.html) function
from the `polycor` package, and assumes factor levels are in the given
order. NAs are automatically handled by pairwise omit.

## Examples

``` r
pair_polychor(iris)
#> Warning: Data has just one column.
```
