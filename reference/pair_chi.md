# Pearson's Contingency Coefficient for association between factors.

Calculates Pearson's Contingency coefficient for every factor variable
pair in a dataset.

## Usage

``` r
pair_chi(d, handle.na = TRUE, warnings = TRUE, ...)
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

A tibble of class `pairwise` with calculated Pearson's contingency
coefficient for every factor variable pair, or NULL if there are not at
least two factor variables

## Details

The Pearson's contingency coefficient is calculated using
[`ContCoef`](https://andrisignorell.github.io/DescTools/reference/CramerV.html).
NAs are automatically handled by pairwise omit.

## Examples

``` r
 pair_chi(iris)
#> Warning: Data has just one column.
```
