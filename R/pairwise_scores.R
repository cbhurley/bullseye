#' Calculates scores or conditional scores for a dataset
#'
#' Calculates scores for every variable pair in a dataset when `by` is `NULL`. If `by`
#' is a name of a variable in the dataset, conditional scores for every
#' variable pair at different levels of the grouping variable are calculated.
#'
#' @param d a dataframe
#' @param by a character string for the name of the conditioning variable. Set to `NULL` by default.
#' @param ungrouped Ignored if `by` is `NULL`.
#'                       If TRUE calculates the ungrouped score in addition to grouped scores.
#' @param control a list for the measures to be calculated for different variable types. The default is
#'              [`pair_control()`] which calculates Pearson's correlation if the variable pair is numeric,
#'              canonical correlation for factor or mixed pairs, and polychoric correlation for two ordered factors.
#' @param handle.na If TRUE uses pairwise complete observations to calculate measure of association.
#' @details Returns a `pairwise` tibble structure.
#' @return A tibble with class `pairwise`.
#'
#' @export
#'
#'
#' @examples
#' irisc <- pairwise_scores(iris)
#' irisc <- pairwise_scores(iris, control=pair_control(nnargs= c(method="spearman")))
#' irisc <- pairwise_scores(iris, control=pair_control(fn="pair_ace"))
#'
#' #Lots of numerical measures
#' irisc <- pairwise_scores(iris, control=pair_control(nn="pairwise_multi", fn=NULL))
#' irisc <- pairwise_scores(iris, 
#'              control=pair_control(nn="pairwise_multi",  nnargs="pair_cor", fn=NULL))

#' #conditional measures
#' cond_iris <- pairwise_scores(iris, by = "Species") 
#' cond_iris_wo <- pairwise_scores(iris, by = "Species",ungrouped=FALSE) # without overall
#' irisc <- pairwise_scores(iris, control=pair_control(nn="pairwise_multi", fn=NULL))
#' irisc <- pairwise_scores(iris, by = "Species",control=pair_control(nn="pairwise_multi", fn=NULL))
#'
#' #scagnostics
#' sc <- pairwise_scores(iris, control=pair_control(nn="pair_scagnostics", fn=NULL)) # ignore fn pairs
#' sc <- pairwise_scores(iris, by = "Species",
#'                   control=pair_control(nn="pair_scagnostics", fn=NULL)) # ignore fn pairs
#' @importFrom rlang .data

pairwise_scores  <- function(d,
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
            if (is.null(entry$argList))
            m <- do.call(get(entry$funName), list(dsub, handle.na=handle.na))
            else m <- do.call(get(entry$funName), list(dsub,  handle.na=handle.na, entry$argList))
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
    measures <- dplyr::bind_rows(measures)
    if (!is.null(control$oo)) {
      ordinals <- which(sapply(d, is.ordered))
      if (length(ordinals) >=2){
        dsub <- d[,ordinals]
        m <- do.call(get(control$oo), c(list(dsub, handle.na=handle.na),control$ooargs))
        measures<-dplyr::rows_update(measures,m,  by = c("x","y"))
      }
    }
    measures
  } else {

    if (!(by %in% names(d))) cli::cli_abort(c("{.var by} not present in dataset."))
    tab <- table(d[[by]])
    if (any(tab == 1)) cli::cli_abort(c("{by} cannot be used as a grouping variable. Need more than one observation at each level."))
    result <- d |>
      dplyr::group_by(.data[[by]]) |>
      dplyr::group_modify(function(x,y) pairwise_scores(x, control=control,handle.na=handle.na)) |>
      dplyr::ungroup() |>
      dplyr::mutate(group=.data[[by]]) |>
      dplyr::select(-dplyr::all_of(by))
    class(result)<-append("pairwise", class(result))
    if (ungrouped){
      overall <- d |>
        dplyr::select(-dplyr::all_of(by)) |>
        pairwise_scores(control=control,handle.na=handle.na)
      result <- rbind(result, overall)
    }
    result

  }

}

#' Default scores calculated  by `pairwise_scores`
#'
#' Gives a list specifying the function to be used for two numeric (nn) variables, two factors (ff), two ordinals (oo)  
#' and for a factor-numeric pair (fn). 
#'
#' @param nn function for numeric pairs of variables, should return object of class `pairwise`. Use NULL to ignore numeric pairs.
#' @param oo function for ordered factor pairs of variables, should return object of class `pairwise`. Use NULL to ignore ordered factor pairs.
#' @param ff function for factor pairs of variables (not ordered), should return object of class `pairwise`. Use NULL to ignore factor-factor pairs.
#' @param fn function for factor-numeric pairs of variables, should return object of class `pairwise`. Use NULL to ignore factor-numeric pairs.
#' @param nnargs other arguments for the nn function
#' @param ooargs other arguments for the oo function
#' @param ffargs other arguments for the ff function
#' @param fnargs other arguments for the fn function
#' 
#' @return list
#' @export
#'

pair_control <- function(nn = "pair_cor",
                         oo = "pair_polychor",
                         ff = "pair_cancor",
                         fn = "pair_cancor",
                         nnargs=NULL,  ooargs=NULL, ffargs=NULL,fnargs=NULL){
  
  list(nn=nn[1], fn=fn[1],oo=oo[1], ff=ff[1],nnargs=nnargs,fnargs=fnargs,ooargs=NULL, ffargs=ffargs)
}
