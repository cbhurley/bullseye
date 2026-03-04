# Goodman Kruskal's Gamma for association between ordinal factors.

Calculates Goodman Kruskal's Gamma coefficient for every factor variable
pair in a dataset.

## Usage

``` r
pair_gkGamma(d, handle.na = TRUE, warnings = TRUE, ...)
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

A tibble of class `pairwise` with factor variable pairs and Goodman
Kruskal's Gamma coefficient, or NULL if there are not at least two
factor variables

## Details

The Goodman Kruskal's Gamma coefficient is calculated using
[`GoodmanKruskalGamma`](https://andrisignorell.github.io/DescTools/reference/GoodmanKruskalGamma.html)
function from the `DescTools` package. Assumes factor levels are in the
given order. NAs are automatically handled by pairwise omit.

## Examples

``` r
 pair_gkGamma(iris)
#> Warning: Data has just one column.
```
