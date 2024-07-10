#' Calculates association or conditional association measures for a dataset
#'
#' Calculates association measures for every variable pair in a dataset when `by` is `NULL`. If `by`
#' is a name of a variable in the dataset, conditional association measures for every
#' variable pair at different levels of the grouping variable are calculated.
#'
#' @param d a dataframe
#' @param by a character string for the name of the conditioning variable. Set to `NULL` by default.
#' @param ungrouped Ignored if `by` is `NULL`.
#'                       If TRUE calculates the overall measure of association for every pair of
#'                       variable, in addition to association measures for pairs at levels of the
#'                       conditioning variable.

#' @param control a list for the measures to be calculated for different variable types. The default is
#'              [`pair_control()`] which calculates Pearson's correlation if the variable pair is numeric,
#'              and canonical correlation for factor or mixed pairs.
#' @param handle.na If TRUE uses pairwise complete observations to calculate measure of association.
#' @details Returns a `pairwise` tibble structure.
#' @return A tibble with class `pairwise`.
#'
#' @export
#'
#'
#' @examples
#' irisc <- pair_scores(iris)
#' irisc <- pair_scores(iris, control=pair_control(nnargs= list(method="spearman")))
#' irisc <- pair_scores(iris, control=pair_control(fn="pair_ace"))
#'
#' #Lots of numerical measures
#' irisc <- pair_scores(iris, control=pair_control(nn="pair_all", fn=NULL))
#' irisc <- pair_scores(iris, control=pair_control(nn=c("pair_all"), nnargs=c("pair_cor"), fn=NULL))

#' #conditional measures
#' cond_iris <- pair_scores(iris, by = "Species") 
#' cond_iris_wo <- pair_scores(iris, by = "Species",ungrouped=FALSE) # without overall
#' irisc <- pair_scores(iris, control=pair_control(nn="pair_all", fn=NULL))
#' irisc <- pair_scores(iris, by = "Species",control=pair_control(nn="pair_all", fn=NULL))
#'
#' #scagnostics
#' sc <- pair_scores(iris, control=pair_control(nn="pair_scagnostics", fn=NULL)) # ignore fn pairs
#' sc <- pair_scores(iris, by = "Species",
#'                   control=pair_control(nn="pair_scagnostics", fn=NULL)) # ignore fn pairs
#' @importFrom rlang .data

pair_scores  <- function(d,
                        by=NULL,
                        ungrouped=TRUE,
                        control=pair_control(),
                        handle.na=TRUE){

  
  check_df(d)
  if(is.null(by)) {
    
   vartypes <- sapply(names(d), function(u)
      if (is.numeric(d[[u]])) "n"
      else if (is.factor(d[[u]])) "f"
      else "other")

    lookup <- function(xy){
      list(funName = control[[xy]],
          argList = control[[paste0(xy,"args")]])
    }
    utypes <- sort(unique(vartypes[vartypes != "other"]))
    
    measures <- vector("list", 3)
    k <- 0
    for (i in seq(along=utypes)){
      for (j in seq(along=utypes)){
        if (i <=j){
          xy <- paste0(utypes[i], utypes[j])
          entry <- lookup(xy)
          dsub <- d[vartypes==utypes[i] | vartypes==utypes[j]]
          if (ncol(dsub)>1 & !is.null(entry$funName)){
            m <- do.call(get(entry$funName), c(list(dsub, handle.na=handle.na), entry$argList))
            if (!inherits(m, "pairwise"))
              stop("Calculated pairwise scores must be of type pairwise")
            m <- m[m$pair_type==xy,]
            k <- k+1
            measures[[k]] <- m[]
            names(measures)[k] <- paste0(utypes[i], utypes[j])
          }
        }
      }
    }
    measures <- measures[1:k]
    dplyr::bind_rows(measures)
  } else {

    if (!(by %in% names(d))) cli::cli_abort(c("{.var by} not present in dataset."))
    tab <- table(d[[by]])
    if (any(tab == 1)) cli::cli_abort(c("{by} cannot be used as a grouping variable. Need more than one observation at each level."))
    # if (!(is.factor(d[[by]]) | is.character(d[[by]])))
    #   cli::cli_abort(c("{.var by} should name a factor or character variable."))
    result <- d |>
      dplyr::group_by(.data[[by]]) |>
      dplyr::group_modify(function(x,y) pair_scores(x, control=control,handle.na=handle.na)) |>
      dplyr::ungroup() |>
      dplyr::mutate(group=.data[[by]]) |>
      dplyr::select(-dplyr::all_of(by))
    class(result)<-append("pairwise", class(result))
    if (ungrouped){
      overall <- d |>
        dplyr::select(-dplyr::all_of(by)) |>
        pair_scores(control=control,handle.na=handle.na)
      result <- rbind(result, overall)
    }
    result

  }

}

#' Default scores calculated  by `pair_scores`
#'
#' Gives a list specifying the function to be used for two numeric (nn) variables, two factors (ff) 
#' and for a factor-numeric pair (fn). 
#'
#' @param nn function for numeric pairs of variables, should return object of class `pairwise`. Use NULL to ignore numeric pairs.
#' @param fn function for factor-numeric pairs of variables, should return object of class `pairwise`. Use NULL to ignore factor-numeric pairs.
#' @param ff function for factor pairs of variables, should return object of class `pairwise`. Use NULL to ignore factor-factor pairs.
#' @param nnargs other arguments for the nn function
#' @param fnargs other arguments for the fn function
#' @param ffargs other arguments for the ff function
#'
#' @return list
#' @export
#'

pair_control <- function(nn = c("pair_cor","pair_dcor","pair_mine","pair_ace",
                                "pair_cancor","pair_nmi", "pair_scagnostics" ),
                         fn = c("pair_cancor", "pair_nmi","pair_ace"),
                         ff = c("pair_cancor","pair_ace","pair_nmi","pair_polycor",    
                                "pair_tau","pair_gkGamma","pair_gkTau","pair_uncertainty",
                                "pair_chi"),
                         nnargs=NULL, fnargs=NULL, ffargs=NULL){
  
  list(nn=nn[1], fn=fn[1],ff=ff[1],nnargs=nnargs,fnargs=fnargs,ffargs=ffargs)
}
