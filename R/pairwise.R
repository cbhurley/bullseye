#' A generic function to create a data structure for every variable pair in a dataset
#'
#' Creates a data structure for every variable pair in a dataset.
#'
#' @param x A dataframe or symmetric matrix.
#' @param score a character string indicating the value of association, either "nn", "fn", "ff".
#' @param pair_type a character string specifying the type of variable pair.
#' @return A tbl_df of class `pairwise` for pairs of variables with a column `value` for the score value,
#' `score` for a type of association value and `pair_type` for the type of variable pair.
#' @export
#'
#' @examples
#' pairwise(cor(iris[,1:4]), score="pearson")
#' pairwise(iris)


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
 
  xindex <- as.vector(row(m))
  yindex <- as.vector(col(m))
  rnames <- rownames(m) %||% paste0("V", seq_along(xindex))
  d <- dplyr::tibble(x=rnames[xindex], y= rnames[yindex],
                     score=score, group="all", value=as.vector(m),pair_type=pair_type)
  class(d)<-append("pairwise", class(d))
  d[d$x<d$y,]
}


#' @describeIn pairwise  pairwise method
#' @export
pairwise.data.frame <- function(x, score=NA_character_, pair_type=NA_character_){
  if(inherits(x, "pairwise")) return(x)
  if (all(identical(names(x),  c("x","y","score","group","value","pair_type"))) && all(x$x < x$y)){
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
                       group="all", value=.data$r, pair_type=NA, .keep="none") |>
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

