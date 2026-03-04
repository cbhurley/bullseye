# Pairwise score functions available in the package

A tibble of score functions along with the types of variable pairs these
functions can be applied to. It also contains information regarding the
packages used to calculate scores and the range of the values
calculated.

## Usage

``` r
pair_methods
```

## Format

An object of class `tbl_df` (inherits from `tbl`, `data.frame`) with 17
rows and 7 columns.

## Value

tibble

## Examples

``` r
 pair_methods
#> # A tibble: 17 × 7
#>    name             nn    ff    fn    from                         range ordinal
#>    <chr>            <lgl> <lgl> <lgl> <chr>                        <chr> <lgl>  
#>  1 pair_cor         TRUE  FALSE FALSE cor                          [-1,… NA     
#>  2 pair_dcor        TRUE  FALSE FALSE energy::dcor2d               [0,1] NA     
#>  3 pair_mine        TRUE  FALSE FALSE minerva::mine                [0,1] NA     
#>  4 pair_ace         TRUE  TRUE  TRUE  acepack::ace                 [0,1] FALSE  
#>  5 pair_cancor      TRUE  TRUE  TRUE  cancor                       [0,1] FALSE  
#>  6 pair_nmi         TRUE  TRUE  TRUE  linkspotter::maxNMI          [0,1] FALSE  
#>  7 pair_polychor    FALSE TRUE  FALSE polycor::polychor            [-1,… TRUE   
#>  8 pair_polyserial  FALSE FALSE TRUE  polycor::polyserial          [-1,… TRUE   
#>  9 pair_tauB        FALSE TRUE  FALSE DescTools::KendalTauB        [-1,… TRUE   
#> 10 pair_tauA        FALSE TRUE  FALSE DescTools::KendalTauA        [-1,… TRUE   
#> 11 pair_tauC        FALSE TRUE  FALSE DescTools::StuartTauC        [-1,… TRUE   
#> 12 pair_tauW        FALSE TRUE  FALSE DescTools::KendalW           [-1,… TRUE   
#> 13 pair_gkGamma     FALSE TRUE  FALSE DescTools::GoodmanKruskalGa… [-1,… TRUE   
#> 14 pair_gkTau       FALSE TRUE  FALSE DescTools::GoodmanKruskalTau [0,1] TRUE   
#> 15 pair_uncertainty FALSE TRUE  FALSE DescTools::UncertCoef        [0,1] FALSE  
#> 16 pair_chi         FALSE TRUE  FALSE DescTools::ContCoef          [0,1] FALSE  
#> 17 pair_scag        TRUE  FALSE FALSE scagnostics::scagnostics     [0,1] NA     
```
