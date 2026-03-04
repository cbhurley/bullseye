# Default scores calculated by `pairwise_scores`

Gives a list specifying the function to be used for two numeric (nn)
variables, two factors (ff), two ordinals (oo) and for a factor-numeric
pair (fn).

## Usage

``` r
pair_control(
  nn = "pair_cor",
  oo = "pair_polychor",
  ff = "pair_cancor",
  fn = "pair_cancor",
  nnargs = NULL,
  ooargs = NULL,
  ffargs = NULL,
  fnargs = NULL
)
```

## Arguments

- nn:

  function for numeric pairs of variables, should return object of class
  `pairwise`. Use NULL to ignore numeric pairs.

- oo:

  function for ordered factor pairs of variables, should return object
  of class `pairwise`. Use NULL to ignore ordered factor pairs.

- ff:

  function for factor pairs of variables (not ordered), should return
  object of class `pairwise`. Use NULL to ignore factor-factor pairs.

- fn:

  function for factor-numeric pairs of variables, should return object
  of class `pairwise`. Use NULL to ignore factor-numeric pairs.

- nnargs:

  other arguments for the nn function

- ooargs:

  other arguments for the oo function

- ffargs:

  other arguments for the ff function

- fnargs:

  other arguments for the fn function

## Value

list
