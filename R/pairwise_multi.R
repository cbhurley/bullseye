

#' Calculates multiple scores
#'
#' Calculates multiple scores for every variable pair in a dataset.
#'
#' @param d dataframe
#' @param scores a character vector naming functions returning a `pairwise` from a dataset.
#'
#' @param handle.na If TRUE uses pairwise complete observations to calculate pairwise score, otherwise NAs not handled.
#' @param warnings If TRUE, generates a warning for datasets of one row, one column, or with constant variables.

#' @return tibble of class "pairwise"
#' @export
#'
#' @examples
#' iris1 <- iris
#' iris1$Sepal.Length <- cut(iris1$Sepal.Length,3)
#' pairwise_multi(iris1)


pairwise_multi <- function(d,scores=c("pair_cor", "pair_spearman", "pair_dcor","pair_mine","pair_ace",
                                    "pair_cancor","pair_nmi","pair_uncertainty",
                                    "pair_chi"),
                         handle.na=TRUE, warnings=TRUE) {
  d <- check_df(d)
  if (!is.character(scores))
    cli::cli_abort(c("{.var scores}  should be of type character."))
  if (warnings) check_constant_var(d)
  scores <- unique(scores)                        
  results <- vector("list", length(scores))
  for (i in 1:length(scores)){
    results[[i]] <- do.call(what = get(scores[i]), args = list(d = d, handle.na = handle.na, warnings=FALSE))
  }

  results <- dplyr::bind_rows(results)

  # if ("pair_cor" %in% scores) {
  #   spearman <- pair_cor(d,method="spearman", handle.na = handle.na, warnings=FALSE)
  #   # kendall <- pair_cor(d,method="kendall", handle.na = handle.na)
  #   results <- dplyr::bind_rows(results, spearman)
  # }

  if ("pair_tau" %in% scores) {
    taua <- pair_tau(d,method="A", handle.na = handle.na, warnings=FALSE)
    tauc <- pair_tau(d,method="C", handle.na = handle.na, warnings=FALSE)
    tauw <- pair_tau(d,method="W", handle.na = handle.na, warnings=FALSE)
    results <- dplyr::bind_rows(results, taua, tauc, tauw)
  }

  results |> 
    dplyr::arrange(.data$x, .data$y)

}
