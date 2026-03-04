# Calculates scores or conditional scores for a dataset

Calculates scores for every variable pair in a dataset when `by` is
`NULL`. If `by` is a name of a variable in the dataset, conditional
scores for every variable pair at different levels of the grouping
variable are calculated.

## Usage

``` r
pairwise_scores(
  d,
  by = NULL,
  ungrouped = TRUE,
  control = pair_control(),
  handle.na = TRUE,
  warnings = TRUE,
  add.nobs = FALSE
)
```

## Arguments

- d:

  a dataframe

- by:

  a character string for the name of the conditioning variable. Set to
  `NULL` by default.

- ungrouped:

  Ignored if `by` is `NULL`. If TRUE calculates the ungrouped score in
  addition to grouped scores.

- control:

  a list for the measures to be calculated for different variable types.
  The default is
  [`pair_control()`](https://cbhurley.github.io/bullseye/reference/pair_control.md)
  which calculates Pearson's correlation if the variable pair is
  numeric, canonical correlation for factor or mixed pairs, and
  polychoric correlation for two ordered factors.

- handle.na:

  If TRUE uses pairwise complete observations to calculate measure of
  association.

- warnings:

  If TRUE, generates a warning for datasets of one row, one column, or
  with constant variables.

- add.nobs:

  If TRUE, adds an extra column containing the number of observations
  used for each score.

## Value

A tibble with class `pairwise`.

## Details

Returns a `pairwise` tibble structure.

## Examples

``` r
irisc <- pairwise_scores(iris)
irisc <- pairwise_scores(iris, control=pair_control(nnargs= c(method="spearman")))
irisc <- pairwise_scores(iris, control=pair_control(fn="pair_ace"))

#Lots of numerical measures
irisc <- pairwise_scores(iris, control=pair_control(nn="pairwise_multi", fn=NULL))
irisc <- pairwise_scores(iris, 
             control=pair_control(nn="pairwise_multi",  nnargs="pair_cor", fn=NULL))
#conditional measures
cond_iris <- pairwise_scores(iris, by = "Species") 
cond_iris_wo <- pairwise_scores(iris, by = "Species",ungrouped=FALSE) # without overall
irisc <- pairwise_scores(iris, control=pair_control(nn="pairwise_multi", fn=NULL))
irisc <- pairwise_scores(iris, by = "Species",control=pair_control(nn="pairwise_multi", fn=NULL))

#scagnostics
sc <- pairwise_scores(iris, control=pair_control(nn="pair_scagnostics", fn=NULL)) # ignore fn pairs
sc <- pairwise_scores(iris, by = "Species",
                  control=pair_control(nn="pair_scagnostics", fn=NULL)) # ignore fn pairs
```
