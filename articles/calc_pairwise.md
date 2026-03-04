# Calculating pairwise scores using bullseye.

`bullseye` is an R package which calculates measures of correlation and
other association scores for pairs of variables in a dataset and offers
visualisations of these measures in different layouts. The package also
calculates and visualises the pairwise scores for different levels of a
grouping variable.

This vignette gives an overview of how these pairwise variable measures
are calculated. Visualisations of these calculated measures are provided
in the accompanying vignette.

``` r
# install.packages("palmerpenguins")

library(bullseye)
library(palmerpenguins)
library(dplyr)
```

Table 1 lists the different measures of association provided in the
package with the variable types they can be used with, the package used
for calculation, the information on whether the measure is symmetric,
and the minimum and maximum value of the measure.

| name             | nn    | ff    | fn    | from                           | range    | ordinal |
|:-----------------|:------|:------|:------|:-------------------------------|:---------|:--------|
| pair_cor         | TRUE  | FALSE | FALSE | cor                            | \[-1,1\] | NA      |
| pair_dcor        | TRUE  | FALSE | FALSE | energy::dcor2d                 | \[0,1\]  | NA      |
| pair_mine        | TRUE  | FALSE | FALSE | minerva::mine                  | \[0,1\]  | NA      |
| pair_ace         | TRUE  | TRUE  | TRUE  | acepack::ace                   | \[0,1\]  | FALSE   |
| pair_cancor      | TRUE  | TRUE  | TRUE  | cancor                         | \[0,1\]  | FALSE   |
| pair_nmi         | TRUE  | TRUE  | TRUE  | linkspotter::maxNMI            | \[0,1\]  | FALSE   |
| pair_polychor    | FALSE | TRUE  | FALSE | polycor::polychor              | \[-1,1\] | TRUE    |
| pair_polyserial  | FALSE | FALSE | TRUE  | polycor::polyserial            | \[-1,1\] | TRUE    |
| pair_tauB        | FALSE | TRUE  | FALSE | DescTools::KendalTauB          | \[-1,1\] | TRUE    |
| pair_tauA        | FALSE | TRUE  | FALSE | DescTools::KendalTauA          | \[-1,1\] | TRUE    |
| pair_tauC        | FALSE | TRUE  | FALSE | DescTools::StuartTauC          | \[-1,1\] | TRUE    |
| pair_tauW        | FALSE | TRUE  | FALSE | DescTools::KendalW             | \[-1,1\] | TRUE    |
| pair_gkGamma     | FALSE | TRUE  | FALSE | DescTools::GoodmanKruskalGamma | \[-1,1\] | TRUE    |
| pair_gkTau       | FALSE | TRUE  | FALSE | DescTools::GoodmanKruskalTau   | \[0,1\]  | TRUE    |
| pair_uncertainty | FALSE | TRUE  | FALSE | DescTools::UncertCoef          | \[0,1\]  | FALSE   |
| pair_chi         | FALSE | TRUE  | FALSE | DescTools::ContCoef            | \[0,1\]  | FALSE   |
| pair_scag        | TRUE  | FALSE | FALSE | scagnostics::scagnostics       | \[0,1\]  | NA      |

List of the functions available in the package for calculating different
association measures along with the packages used for calculation.

## Calculating correlation and other association measures

Each of the functions in the first column of Table 1 calculates pairwise
scores for a dataset.

``` r
sc_dcor <- pair_dcor(penguins)
str(sc_dcor)
#> pairwise [10 × 6] (S3: pairwise/tbl_df/tbl/data.frame)
#>  $ x        : chr [1:10] "bill_depth_mm" "bill_length_mm" "bill_depth_mm" "body_mass_g" ...
#>  $ y        : chr [1:10] "bill_length_mm" "flipper_length_mm" "flipper_length_mm" "flipper_length_mm" ...
#>  $ score    : chr [1:10] "dcor" "dcor" "dcor" "dcor" ...
#>  $ group    : chr [1:10] "all" "all" "all" "all" ...
#>  $ value    : Named num [1:10] 0.387 0.666 0.704 0.867 0.587 ...
#>   ..- attr(*, "names")= chr [1:10] "" "" "" "" ...
#>  $ pair_type: chr [1:10] "nn" "nn" "nn" "nn" ...
```

For example, we see that `pair_dcor` calculates the distance correlation
for every pair of numeric variables in the `penguins` dataset. There are
missing values in this dataset, all the `pair_` functions use pairwise
complete observations by default.

`sc_dcor` is a tibble of class `pairwise`, with the two variables in
columns `x` and `y` (arranged in alphabetical order), calculated values
in the column `value`, and the name of the score calculated in the
column `score`. All of the variables are numeric, hence “nn” in the
`pair_type` column.

Similarly, one can use `pair_nmi` to calculate normalised mutual
information for numeric, factor and mixed pairs of variables.

``` r
sc_nmi <- pair_nmi(penguins)
sc_nmi
#> # A tibble: 28 × 6
#>    x                 y                 score group  value pair_type
#>    <chr>             <chr>             <chr> <chr>  <dbl> <chr>    
#>  1 bill_depth_mm     bill_length_mm    nmi   all   0.225  nn       
#>  2 bill_length_mm    flipper_length_mm nmi   all   0.375  nn       
#>  3 bill_depth_mm     flipper_length_mm nmi   all   0.470  nn       
#>  4 body_mass_g       flipper_length_mm nmi   all   0.581  nn       
#>  5 bill_length_mm    body_mass_g       nmi   all   0.303  nn       
#>  6 bill_depth_mm     body_mass_g       nmi   all   0.443  nn       
#>  7 bill_length_mm    year              nmi   all   0.0517 nn       
#>  8 bill_depth_mm     year              nmi   all   0.0387 nn       
#>  9 flipper_length_mm year              nmi   all   0.0707 nn       
#> 10 body_mass_g       year              nmi   all   0.0445 nn       
#> # ℹ 18 more rows
```

The main difference here is that factor variables are included. In the
`pair_type` column, “ff” and “fn” indicate factor-factor and
factor-numeric pairs.

If you want more control over the measure calculated, the function
`pairwise_scores` calculates a different score depending on variable
types.

``` r
pairwise_scores(penguins) |> distinct(score, pair_type)
#> # A tibble: 3 × 2
#>   score   pair_type
#>   <chr>   <chr>    
#> 1 cancor  ff       
#> 2 cancor  fn       
#> 3 pearson nn
```

As you can see, the default uses pearson’s correlation for numeric
pairs, and canonical correlation for factor-numeric or factor-factor
pairs. In addition polychoric correlation is used for two ordered
factors, but there are no ordered factors in this data. Alternative
scores may be specified using the `control` argument to
`pairwise_scores`. The default value for this `control` argument is
given by

``` r
pair_control()
#> $nn
#> [1] "pair_cor"
#> 
#> $fn
#> [1] "pair_cancor"
#> 
#> $oo
#> [1] "pair_polychor"
#> 
#> $ff
#> [1] "pair_cancor"
#> 
#> $nnargs
#> NULL
#> 
#> $fnargs
#> NULL
#> 
#> $ooargs
#> NULL
#> 
#> $ffargs
#> NULL
```

## Calculating multiple measures

If you want for instance to compare distance correlation and mutual
information measures in a display, two `pairwise` data structures can be
combined:

``` r
bind_rows(sc_dcor, sc_nmi) |> arrange(x,y)
#> # A tibble: 38 × 6
#>    x             y                 score group value pair_type
#>    <chr>         <chr>             <chr> <chr> <dbl> <chr>    
#>  1 bill_depth_mm bill_length_mm    dcor  all   0.387 nn       
#>  2 bill_depth_mm bill_length_mm    nmi   all   0.225 nn       
#>  3 bill_depth_mm body_mass_g       dcor  all   0.614 nn       
#>  4 bill_depth_mm body_mass_g       nmi   all   0.443 nn       
#>  5 bill_depth_mm flipper_length_mm dcor  all   0.704 nn       
#>  6 bill_depth_mm flipper_length_mm nmi   all   0.470 nn       
#>  7 bill_depth_mm island            nmi   all   0.282 fn       
#>  8 bill_depth_mm sex               nmi   all   0.356 fn       
#>  9 bill_depth_mm species           nmi   all   0.493 fn       
#> 10 bill_depth_mm year              dcor  all   0.112 nn       
#> # ℹ 28 more rows
```

We provide another function `pairwise_multi` which calculates multiple
association measures for every variable pair in a dataset. By default
this function combines the results of `pair_cor`,
`pair_dcor`,`pair_mine`,`pair_ace`,
`pair_cancor`,`pair_nmi`,`pair_uncertainty`, `pair_chi`, but any subset
of the `pair_` functions may be supplied as an argument, as in the
second example below.

``` r
pairwise_multi(penguins)
#> # A tibble: 130 × 6
#>    x             y              score    group  value pair_type
#>    <chr>         <chr>          <chr>    <chr>  <dbl> <chr>    
#>  1 bill_depth_mm bill_length_mm pearson  all   -0.235 nn       
#>  2 bill_depth_mm bill_length_mm spearman all   -0.222 nn       
#>  3 bill_depth_mm bill_length_mm dcor     all    0.387 nn       
#>  4 bill_depth_mm bill_length_mm MIC      all    0.313 nn       
#>  5 bill_depth_mm bill_length_mm ace      all    0.585 nn       
#>  6 bill_depth_mm bill_length_mm cancor   all    0.235 nn       
#>  7 bill_depth_mm bill_length_mm nmi      all    0.225 nn       
#>  8 bill_depth_mm body_mass_g    pearson  all   -0.472 nn       
#>  9 bill_depth_mm body_mass_g    spearman all   -0.432 nn       
#> 10 bill_depth_mm body_mass_g    dcor     all    0.614 nn       
#> # ℹ 120 more rows
dcor_nmi <- pairwise_multi(penguins, c("pair_dcor", "pair_nmi"))
```

## Calculating grouped measures

For each of the pairwise calculation functions, they can be wrapped
using `pairwise_by` to build a score calculation for each level of a
grouping variable. Of course, grouped scores could be calculated using
`dplyr` machinery, but it is a bit more work.

``` r
pairwise_by(penguins, by="species", pair_cor)
#> Warning: Variable island of data for group species = Chinstrap has at most one unique
#> value, discounting NA.
#> Warning: Variable island of data for group species = Gentoo has at most one unique
#> value, discounting NA.
#> # A tibble: 40 × 6
#>    x             y                 score   group      value pair_type
#>    <chr>         <chr>             <chr>   <fct>      <dbl> <chr>    
#>  1 bill_depth_mm bill_length_mm    pearson Adelie     0.391 nn       
#>  2 bill_depth_mm bill_length_mm    pearson Chinstrap  0.654 nn       
#>  3 bill_depth_mm bill_length_mm    pearson Gentoo     0.643 nn       
#>  4 bill_depth_mm bill_length_mm    pearson all       -0.235 nn       
#>  5 bill_depth_mm body_mass_g       pearson Adelie     0.576 nn       
#>  6 bill_depth_mm body_mass_g       pearson Chinstrap  0.604 nn       
#>  7 bill_depth_mm body_mass_g       pearson Gentoo     0.719 nn       
#>  8 bill_depth_mm body_mass_g       pearson all       -0.472 nn       
#>  9 bill_depth_mm flipper_length_mm pearson Adelie     0.308 nn       
#> 10 bill_depth_mm flipper_length_mm pearson Chinstrap  0.580 nn       
#> # ℹ 30 more rows
```

Use argument `ungrouped=FALSE` to suppress calculation of the ungrouped
scores.

The `pairwise_scores` function introduced previously also has a `by`
argument, and provides pairwise scores for the levels of a grouping
variable.

``` r
sc_sex <- pairwise_scores(penguins, by="species")
#> Warning: Variable island of data for group species = Chinstrap has at most one unique
#> value, discounting NA.
#> Warning: Variable island of data for group species = Gentoo has at most one unique
#> value, discounting NA.
```

The column `group` now shows the levels of the grouping variable, along
with “all” for ungrouped scores. Use `ungrouped=FALSE` to suppress
calculation of the ungrouped scores.

``` r
sc_sex |> distinct(group)
#> # A tibble: 4 × 1
#>   group    
#>   <fct>    
#> 1 Adelie   
#> 2 Chinstrap
#> 3 Gentoo   
#> 4 all
```

If you want to calculate different scores to the default, specify this
via the `control` argument:

``` r
pc <- pair_control(nn="pairwise_multi", nnargs= c("pair_dcor", "pair_ace"), fn=NULL, ff=NULL)
sc_sex <- pairwise_scores(penguins, by="species", control=pc, ungrouped=FALSE) 
#> Warning: Variable island of data for group species = Chinstrap has at most one unique
#> value, discounting NA.
#> Warning: Variable island of data for group species = Gentoo has at most one unique
#> value, discounting NA.
```

Both of the functions `pairwise_by` and `pairwise_scores` have an
additional argument `add.nobs`. When this is set to TRUE, the pairwise
structure has an additional column giving the number of observations
used in the calculation of each score.

This is potentially useful information in the presence of imbalanced
data or with missing values.

``` r
pairwise_scores(penguins, by="species",add.nobs=TRUE) 
#> Warning: Variable island of data for group species = Chinstrap has at most one unique
#> value, discounting NA.
#> Warning: Variable island of data for group species = Gentoo has at most one unique
#> value, discounting NA.
#> # A tibble: 84 × 7
#>    x                 y      score  group   value pair_type     n
#>    <chr>             <chr>  <chr>  <fct>   <dbl> <chr>     <int>
#>  1 island            sex    cancor Adelie 0.0164 ff          146
#>  2 bill_length_mm    island cancor Adelie 0.0838 fn          151
#>  3 bill_depth_mm     island cancor Adelie 0.0629 fn          151
#>  4 flipper_length_mm island cancor Adelie 0.148  fn          151
#>  5 body_mass_g       island cancor Adelie 0.0208 fn          151
#>  6 bill_length_mm    sex    cancor Adelie 0.590  fn          146
#>  7 bill_depth_mm     sex    cancor Adelie 0.597  fn          146
#>  8 flipper_length_mm sex    cancor Adelie 0.355  fn          146
#>  9 body_mass_g       sex    cancor Adelie 0.738  fn          146
#> 10 island            year   cancor Adelie 0.104  fn          152
#> # ℹ 74 more rows
```

For other `pairwise` constructors, the `n` column can be added by
calling the `add_nobs_to_pairwise` function directly:

``` r
pair_cancor(penguins) |>
  add_nobs_to_pairwise(penguins) |> pull(n)
#>  [1] 344 342 342 342 342 333 342 342 342 342 342 342 342 342 342 342 333 333 333
#> [20] 333 333 344 344 342 342 342 342 333
```

## Scagnostic measures

The package `scagnostics` provides pairwise variable scores based on
graph-theoretic interestingness measures, for numeric variable pairs
only.

``` r
pair_scagnostics(penguins[,1:5], scagnostic=c("Stringy", "Clumpy"))
#> # A tibble: 6 × 6
#>   x              y                 score   group  value pair_type
#>   <chr>          <chr>             <chr>   <chr>  <dbl> <chr>    
#> 1 bill_depth_mm  bill_length_mm    Stringy all   0.331  nn       
#> 2 bill_depth_mm  bill_length_mm    Clumpy  all   0.0328 nn       
#> 3 bill_depth_mm  flipper_length_mm Stringy all   0.378  nn       
#> 4 bill_depth_mm  flipper_length_mm Clumpy  all   0.530  nn       
#> 5 bill_length_mm flipper_length_mm Stringy all   0.370  nn       
#> 6 bill_length_mm flipper_length_mm Clumpy  all   0.0388 nn
```

Note that the first two variables of the penguins data are non-numeric
and so are ignored in the above calculation.

For group-wise calculation:

``` r
pairwise_by(penguins[,1:5], by="species",function(x) pair_scagnostics(x, scagnostic=c("Stringy", "Clumpy")))
#> Warning: Variable island of data for group species = Chinstrap has at most one unique
#> value, discounting NA.
#> Warning: Variable island of data for group species = Gentoo has at most one unique
#> value, discounting NA.
#> # A tibble: 24 × 6
#>    x             y                 score   group      value pair_type
#>    <chr>         <chr>             <chr>   <fct>      <dbl> <chr>    
#>  1 bill_depth_mm bill_length_mm    Stringy Adelie    0.278  nn       
#>  2 bill_depth_mm bill_length_mm    Clumpy  Adelie    0.0477 nn       
#>  3 bill_depth_mm bill_length_mm    Stringy Chinstrap 0.325  nn       
#>  4 bill_depth_mm bill_length_mm    Clumpy  Chinstrap 0.0758 nn       
#>  5 bill_depth_mm bill_length_mm    Stringy Gentoo    0.393  nn       
#>  6 bill_depth_mm bill_length_mm    Clumpy  Gentoo    0.0579 nn       
#>  7 bill_depth_mm bill_length_mm    Stringy all       0.331  nn       
#>  8 bill_depth_mm bill_length_mm    Clumpy  all       0.0328 nn       
#>  9 bill_depth_mm flipper_length_mm Stringy Adelie    0.392  nn       
#> 10 bill_depth_mm flipper_length_mm Clumpy  Adelie    0.0598 nn       
#> # ℹ 14 more rows
```

## Converting symmetric matrices to `pairwise` and vice-versa.

The conventional way of representing pairwise scores or correlations is
via a numeric symmetric matrix. The tidy `pairwise` structure we use in
`bullseye` is more flexible, and is amenable to multiple measures of
association and grouped measures.

It is straightforward to convert from a symmetric matrix to `pairwise`:

``` r
x <- cor(penguins[, c("bill_length_mm", "bill_depth_mm" ,"flipper_length_mm" ,"body_mass_g")], 
         use= "pairwise.complete.obs")
pairwise(x, score="pearson", pair_type = "nn")
#> # A tibble: 6 × 6
#>   x              y                 score   group  value pair_type
#>   <chr>          <chr>             <chr>   <chr>  <dbl> <chr>    
#> 1 bill_depth_mm  bill_length_mm    pearson all   -0.235 nn       
#> 2 bill_length_mm flipper_length_mm pearson all    0.656 nn       
#> 3 bill_depth_mm  flipper_length_mm pearson all   -0.584 nn       
#> 4 body_mass_g    flipper_length_mm pearson all    0.871 nn       
#> 5 bill_length_mm body_mass_g       pearson all    0.595 nn       
#> 6 bill_depth_mm  body_mass_g       pearson all   -0.472 nn
```

And for the reverse, converting a `pairwise` to a symmetric matrix:

``` r
as.matrix(sc_dcor)
#>                   bill_depth_mm bill_length_mm body_mass_g flipper_length_mm
#> bill_depth_mm                NA     0.38720211   0.6141631         0.7039636
#> bill_length_mm        0.3872021             NA   0.5871319         0.6664558
#> body_mass_g           0.6141631     0.58713186          NA         0.8674122
#> flipper_length_mm     0.7039636     0.66645577   0.8674122                NA
#> year                  0.1117057     0.07842516   0.0790560         0.1643876
#>                         year
#> bill_depth_mm     0.11170568
#> bill_length_mm    0.07842516
#> body_mass_g       0.07905600
#> flipper_length_mm 0.16438763
#> year                      NA
```

### Converting structures from package `correlation`:

`correlation` package calculates different kinds of correlations, such
as partial correlations, Bayesian correlations, multilevel correlations,
polychoric correlations, biweight, percentage bend or Sheperd’s Pi
correlations, distance correlation and more. The output data structure
is a tidy dataframe with a correlation value and correlation tests for
variable pairs for which the correlation method is defined.

``` r
correlation::correlation(penguins)
#> # Correlation Matrix (pearson-method)
#> 
#> Parameter1        |        Parameter2 |     r |         95% CI | t(340) |         p
#> -----------------------------------------------------------------------------------
#> bill_length_mm    |     bill_depth_mm | -0.24 | [-0.33, -0.13] |  -4.46 | < .001***
#> bill_length_mm    | flipper_length_mm |  0.66 | [ 0.59,  0.71] |  16.03 | < .001***
#> bill_length_mm    |       body_mass_g |  0.60 | [ 0.52,  0.66] |  13.65 | < .001***
#> bill_length_mm    |              year |  0.05 | [-0.05,  0.16] |   1.01 | 0.797    
#> bill_depth_mm     | flipper_length_mm | -0.58 | [-0.65, -0.51] | -13.26 | < .001***
#> bill_depth_mm     |       body_mass_g | -0.47 | [-0.55, -0.39] |  -9.87 | < .001***
#> bill_depth_mm     |              year | -0.06 | [-0.17,  0.05] |  -1.11 | 0.797    
#> flipper_length_mm |       body_mass_g |  0.87 | [ 0.84,  0.89] |  32.72 | < .001***
#> flipper_length_mm |              year |  0.17 | [ 0.06,  0.27] |   3.17 | 0.007**  
#> body_mass_g       |              year |  0.04 | [-0.06,  0.15] |   0.78 | 0.797    
#> 
#> p-value adjustment method: Holm (1979)
#> Observations: 342
```

The default calculation uses Pearson correlation. Other options are
available via the `method` argument.

As there is an `as.matrix` method provided for the results of
[`correlation::correlation`](https://easystats.github.io/correlation/reference/correlation.html),
it is straightforward to convert this to a `pairwise` tibble.

``` r
x <- correlation::correlation(penguins)
pairwise(as.matrix(x)) 
#> # A tibble: 10 × 6
#>    x                 y                 score group   value pair_type
#>    <chr>             <chr>             <chr> <chr>   <dbl> <chr>    
#>  1 bill_depth_mm     bill_length_mm    NA    all   -0.235  NA       
#>  2 bill_length_mm    flipper_length_mm NA    all    0.656  NA       
#>  3 bill_depth_mm     flipper_length_mm NA    all   -0.584  NA       
#>  4 body_mass_g       flipper_length_mm NA    all    0.871  NA       
#>  5 bill_length_mm    body_mass_g       NA    all    0.595  NA       
#>  6 bill_depth_mm     body_mass_g       NA    all   -0.472  NA       
#>  7 bill_length_mm    year              NA    all    0.0545 NA       
#>  8 bill_depth_mm     year              NA    all   -0.0604 NA       
#>  9 flipper_length_mm year              NA    all    0.170  NA       
#> 10 body_mass_g       year              NA    all    0.0422 NA
```
