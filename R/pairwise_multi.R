

#' Calculates multiple scores
#'
#' Calculates multiple scores for every variable pair in a dataset.
#'
#' @param d dataframe
#' @param scores a vector naming functions returning a `pairwise` from a dataset.
#'
#' @param handle.na If TRUE uses pairwise complete observations to calculate pairwise score, otherwise NAs not handled.
#'
#' @return tibble of class "pairwise"
#' @export
#'
#' @examples
#' iris1 <- iris
#' iris1$Sepal.Length <- cut(iris1$Sepal.Length,3)
#' pairwise_multi(iris1)


pairwise_multi <- function(d,scores=c("pair_cor", "pair_dcor","pair_mine","pair_ace",
                                    "pair_cancor","pair_nmi","pair_uncertainty",
                                    "pair_chi"),
                         handle.na=TRUE) {
  check_df(d)
  results <- vector("list", length(scores))
  for (i in 1:length(scores)){
    results[[i]] <- do.call(what = get(scores[i]), args = list(d = d, handle.na = handle.na))
  }

  results <- dplyr::bind_rows(results)

  if ("pair_cor" %in% scores) {
    spearman <- pair_cor(d,method="spearman", handle.na = handle.na)
    # kendall <- pair_cor(d,method="kendall", handle.na = handle.na)
    results <- dplyr::bind_rows(results, spearman)
  }

  if ("pair_tau" %in% scores) {
    taua <- pair_tau(d,method="A", handle.na = handle.na)
    tauc <- pair_tau(d,method="C", handle.na = handle.na)
    tauw <- pair_tau(d,method="W", handle.na = handle.na)
    results <- dplyr::bind_rows(results, taua, tauc, tauw)
  }

  results |> 
    dplyr::arrange(.data$x, .data$y)

}
