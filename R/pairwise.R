#' A generic function to create a data structure for summarising  variable pairs in a dataset
#'
#' Creates a data structure for every variable pair in a dataset.
#'
#' @param x A dataframe or symmetric matrix.
#' @param score a character string indicating the value of association
#' @param pair_type a character string specifying the type of variable pair, should be either "nn", "fn", "ff", 
#' for a numeric-numeric pair, factor-numeric pair, or factor-factor pair, or NA if unknown.
#' @return A tbl_df of class `pairwise` for pairs of variables with a column `value` for the score value,
#' `score` for a type of association value and `pair_type` for the type of variable pair.
#' @details The `pairwise` class has columns x and y for (ordered pairs) of variables, where x < y. 
#' The column score has the name of the summary measure used for the two variables, 
#' and the column value has the associated value.
#' The group column defaults to "all", meaning summary measures apply to the complete dataset, 
#' otherwise it describes a subset of the data.
#' The functions `pair_*` calculate pairwise tibbles for the summary measure named by `*`, eg [pair_cor()], [pair_cancor()].
#' The functions [pairwise_scores()] and [pairwise_by()] calculate pairwise tibbles for levels of a grouping variable.
#' The function [pairwise_multi()] calculates a pairwise_tibble for multiple named scores. 
#' The pairwise tibble has at most one row for each combination of x, y, score and group.
#' This is checked prior to plotting by `plot.pairwise`.
#' Note that the pair_type column is included for information purposes, but it is not currently used by `plot.pairwise`.
#' @export
#'
#' @examples
#' pairwise(cor(iris[,1:4]), score="pearson")
#' pairwise(iris)
#' pair_cor(iris)
#' pair_cancor(iris)
#' pairwise_scores(iris, by="Species")
#' pairwise_multi(iris)

pairwise <- function(x, score=NA_character_, pair_type=NA_character_){
  UseMethod("pairwise", x)
}


pairwise_to_matrix <- function(scores, stat=dplyr::first, default=NA){
  allvars <- unique(c(scores$x, scores$y))
  
  scores1 <- dplyr::summarise(scores, 
                               measure= stat(.data$value),
                              .by=dplyr::all_of(c("x","y")))
  scores1 <- scores1[!is.na(scores1$measure),]
  m <- matrix(default, nrow=length(allvars), ncol=length(allvars))
  rownames(m)<- colnames(m)<- allvars
  m[cbind(scores1$x,scores1$y)]<- m[cbind(scores1$y,scores1$x)]<-scores1$measure
  m
}

#' @describeIn pairwise  pairwise method
#' @export
pairwise.matrix <- function(x, score=NA_character_, pair_type=NA_character_){
  m <- x
  if (!isSymmetric(m))
    stop("Input must be a symmetric matrix")
  if (!(pair_type %in% c(NA, "nn", "ff", "fn")))
    cli::cli_warn(c('{.var pair_type} should be in c(NA, "nn", "ff", "fn") '))
  xindex <- as.vector(row(m))
  yindex <- as.vector(col(m))
  rnames <- rownames(m) 
  if (is.null(rnames)) rnames <-  paste0("V", seq_along(xindex))
  d <- dplyr::tibble(x=rnames[xindex], y= rnames[yindex],
                     score=score, group="all", value=as.vector(m),pair_type=pair_type)
  class(d)<-append("pairwise", class(d))
  d[d$x<d$y,]
}


#' @describeIn pairwise  pairwise method
#' @export
pairwise.data.frame <- function(x, score=NA_character_, pair_type=NA_character_){
  if(inherits(x, "pairwise")) return(x)
  if (all(identical(names(x)[1:6],  c("x","y","score","group","value","pair_type"))) && all(x$x < x$y)){
    maxg <- x |> dplyr::summarise(n=dplyr::n(), .by=c("x","y","score","group")) |> dplyr::pull(.data$n) |> max()
    if (maxg==1) {
      class(x)<- c("pairwise", class(x))
     return(x)
    }
  }
  dcor <- diag(ncol(x)) 
  dcor[]<- NA
  rownames(dcor)<- colnames(dcor) <- names(x)
  dcor <- pairwise.matrix(dcor, score=score, pair_type=pair_type)
  if (is.na(pair_type)){
    fn_pair_type <- function(u,v){
      if(is.numeric(x[[u]]) & is.numeric(x[[v]])) {
        "nn"
      } else if(is.factor(x[[u]]) & is.factor(x[[v]])){
        "ff"
      } else {
        "fn"
      }
    }
    dcor$pair_type <- mapply(fn_pair_type, dcor$x,dcor$y, USE.NAMES = FALSE)
 
  }
  dcor
}

#' @describeIn pairwise  pairwise method
#' @export
pairwise.easycorrelation <- function(x, score=NA_character_, pair_type=NA_character_){
  res <- dplyr::mutate(x, x=pmin(.data$Parameter1, .data$Parameter2),
                       y=pmax(.data$Parameter1, .data$Parameter2),
                       score=.data$Method, 
                       group="all", value=.data$r, pair_type=pair_type, .keep="none") |>
    dplyr::filter(.data$x != .data$y)
  class(res) <- class(as.pairwise(diag(1)))
  res <- unique(res)
  res
}

#' @describeIn pairwise  Same as `pairwise`
#' @export
as.pairwise <- function(x, score=NA_character_, pair_type=NA_character_){
  pairwise(x, score=score, pair_type=pair_type)
}

#' Converts a pairwise to a symmetric matrix. Uses the first entry for each (x,y) pair.
#'
#' @param x An object of class pairwise
#' @param ... other arguments
#' @return A symmetric matrix
#' @export
#' 
as.matrix.pairwise <- function(x, ...){
  pairwise_to_matrix(x, stat=dplyr::first,...)
}



#' Adds number of observations column to pairwise tibble
#'
#' @param scores An object of class `pairwise`, calculated from dataset d.
#' @param d a dataframe
#' @param by a character string for the name of the conditioning variable from d used to construct scores. Set to `NULL` by default.
#'
#' @returns A pairwise tibble with a column `n`
#' @export
#'
#' @examples
#' irisc <- pairwise_scores(iris[40:150,], by= "Species")
#' irisc <- add_nobs_to_pairwise(irisc, iris[40:150,], by= "Species")
#' irisc
#' plot_pairwise(irisc) # setosa gets a small slice in proportion to n
add_nobs_to_pairwise <- function(scores, d, by=NULL){
  check_pairwise(scores)
  if (is.null(by))
    scores_new <- 
      scores |> 
      dplyr::rowwise() |>
      dplyr::mutate(n = nrow(na.omit(d[, c(.data$x, .data$y)]))) |>
      dplyr::ungroup()
  else {
    subsets <- split(d,d[[by]], drop=TRUE)
    scores_new <- 
      scores |>
      dplyr::rowwise() |>
      dplyr::mutate(n = ifelse(.data[["group"]] =="all", 
                           nrow(na.omit(d[, c(.data$x, .data$y)])),
                           nrow(na.omit(subsets[[.data[["group"]]]][, c(.data$x, .data$y)])))) |>
      dplyr::ungroup()
  }
   class(scores_new) <- class(scores)
   scores_new
}