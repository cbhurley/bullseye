# Calculates ace based transformations and correlation, handling missing values and factors.

Calculates ace based transformations and correlation, handling missing
values and factors.

## Usage

``` r
ace_cor(x, y, handle.na = TRUE)
```

## Arguments

- x:

  a numeric vector or factor

- y:

  a numeric vector or factor

- handle.na:

  If TRUE uses pairwise complete observations.

## Value

result of acepack::ace

## Examples

``` r
ace_cor(iris$Sepal.Length, iris$Species)
#> 
#> Alternating Conditional Expections
#> 
#> p = 1 , N = 150 
#> 
#> Raw Multiple R-squared: 0.6124 
#> Transformed Multiple R-squared: 0.7028 
#> 
```
