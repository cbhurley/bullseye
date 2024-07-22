
check_df <- function(d){
  if(!inherits(d, "data.frame")){
    cli::cli_abort(c("{.var d} should be a dataframe or tibble."))
  }
}


check_pairwise <- function(d){
   if(!inherits(d, "pairwise")){
     if (!identical(names(d), names(pairwise(data.frame()))))
     cli::cli_abort(c("{.var d} should be of class pairwise."))
  }
}