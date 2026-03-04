# Plot method for class `pairwise`.

Plot method for class `pairwise`.

## Usage

``` r
# S3 method for class 'pairwise'
plot(x, type = c("matrix", "linear"), ...)
```

## Arguments

- x:

  An object of class `pairwise`

- type:

  If "matrix", calls `plot_pairwise`, if "linear" calls
  `plot_pairwise_linear`

- ...:

  further arguments to `plot_pairwise` or `plot_pairwise_linear`

## Value

a plot

## Examples

``` r
plot(pairwise_scores(iris))
```
