
check_df <- function(d){
  if(!inherits(d, "data.frame")){
    cli::cli_abort(c("{.var d} should be a dataframe or tibble."))
  }
}
