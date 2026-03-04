# Uncertainty coefficient for association between factors.

Calculates uncertainty coefficient for every factor variable pair in a
dataset.

## Usage

``` r
pair_uncertainty(d, handle.na = TRUE, warnings = TRUE, ...)
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

A tibble of class `pairwise` with every factor variable pair and
uncertainty coefficient value, or NULL if there are not at least two
factor variables

## Details

The Uncertainty coefficient is calculated using
[`UncertCoef`](https://andrisignorell.github.io/DescTools/reference/UncertCoef.html)
function from the `DescTools` package.

## Examples

``` r
 pair_uncertainty(iris)
#> Warning: Data has just one column.
```
