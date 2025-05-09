
check_df <- function(d){
  if(!inherits(d, "data.frame")){
    cli::cli_abort(c("{.var d} should be a dataframe or tibble."))
  }
  if(inherits(d, "grouped_df")){
    cli::cli_warn(c("Grouping structure in {.var d} is ignored."))
    class(d) <- class(d)[class(d) != "grouped_df"]
  }
  d
}

check_constant_var <- function(d, msg=NULL){
  if (nrow(d) ==1) 
    cli::cli_warn(c("Data has just one row"))
  if (ncol(d) ==1) 
    cli::cli_warn(c("Data has just one column."))
  else for (i in seq(along=d)){
    if (length(na.omit(unique(d[[i]]))) <=1){
      if (is.null(msg)) msg <- "Variable {names(d)[i]} of data has at most one unique value, discounting NA."
      cli::cli_warn(c(msg))
    }
  }
}

# 
# check_pairwise <- function(d){
#    if(!inherits(d, "pairwise")){
#      if (!identical(names(d), names(pairwise(data.frame()))))
#      cli::cli_abort(c("{.var d} should be of class pairwise."))
#   }
# }


check_pairwise <- function(d){
  if(!inherits(d, "pairwise")){
    nm <- c("x","y","score","group","value","pair_type")
    if (! all(nm %in% names(d)))
      cli::cli_abort(paste0("{.var d} must have columns named ", paste0(nm, collapse=","),collapse=""))
  }
  if (!is.character(d$x))
    cli::cli_abort(c("Column {.var x} of  {.var d} should be of type character."))
  if (!is.character(d$y))
    cli::cli_abort(c("Column {.var y} of  {.var d} should be of type character."))
  if (!all(d$x<d$y))
    cli::cli_abort(c("Values in column {.var x} should be less than those in column {.var y}."))
  if (!is.character(d$group) & !is.factor(d$group))
    cli::cli_abort(c("Column {.var group} of  {.var d} should be of type character or a factor."))
  if (!is.character(d$score))
    cli::cli_abort(c("Column {.var score} of  {.var d} should be of type character."))
  if (!is.numeric(d$value))
    cli::cli_abort(c("Column {.var value} of  {.var d} should be of type numeric."))
  if (!is.character(d$pair_type))
    cli::cli_abort(c("Column {.var pair_type} of  {.var d} should be of type character."))
  maxg <- d |> 
    dplyr::summarise(n=dplyr::n(), .by=c("x","y","score","group")) |> 
    dplyr::pull(.data$n) |> max()
  if (maxg > 1) cli::cli_abort(c("{.var d} should have at most one row for each unique combination of 
                                 {.var x}, {.var y}, {.var score},{.var group}"))
  return(TRUE)
}


subset_checks <- function(d, by, warnings=TRUE){
  if (!isTRUE(length(by) == 1)) cli::cli_abort(c("{.var by} must have length 1."))
  if (!(by %in% names(d))) cli::cli_abort(c("{.var by} not present in dataset."))
  tab <- table(d[[by]])
  if (any(tab == 1)) cli::cli_abort(c("{by} cannot be used as a grouping variable. Need more than one observation at each level."))
  if (warnings) {
    subsets <- split(d,d[[by]], drop=TRUE)
    
    for (j in seq(along=subsets)){
      d <- subsets[[j]]
       d[[by]] <- NULL
      for (i in seq(along=d)){
        if (length(na.omit(unique(d[[i]]))) <=1){
          msg <- paste0("Variable {names(d)[i]} of data for group ", by, " = ",names(subsets)[j], " has at most one unique value, discounting NA.")
          cli::cli_warn(c(msg))
        }
      }
    }
  }
}



