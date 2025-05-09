# bullseye 1.0.0

* added functions pair_spearman, pair_kendall
* added pal argument to plot_pairwise and plot_pairwise_linear
* added add.nobs argument to pairwise_by and pairwise_scores
* improved documentation for function pairwise
* new function add_nobs_to_pairwise
* Added warnings to pairwise calculations for variables that are constant.
* Added checks that pairwise object has the required structure prior to plotting

# bullseye 0.1.2

* Minor change to pair_tauW for compatibility with DescTools::KendallW version 0.99.59 handling of NAs.
* Re-instated pair_ace test
* Added dependency on R >= 4.1.0


# bullseye 0.1.1

* Minor fix to pairwise.matrix
* Temporarily remove test using acepack until that package is updated


# bullseye 0.1.0

* Initial CRAN submission.
