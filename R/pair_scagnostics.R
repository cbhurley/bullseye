#' Graph-theoretic scagnostics values
#'
#' Calculates  scagnostic values for every numeric variable pair in a dataset.
#'
#' @param d A dataframe
#' @param scagnostic a character vector for the scagnostic to be calculated. Subset of "Outlying",
#' "Stringy", "Striated", "Clumpy", "Sparse", "Skewed", "Convex", "Skinny" or "Monotonic"
#' @param handle.na If TRUE uses pairwise complete observations.
#' @param ... other arguments
#'
#' @return A tibble of class `pairwise` with  scagnostic values for every numeric variable pair,
#'  or NULL if there are not at least two numeric variables

#' @details The scagnostic values are calculated using \code{\link[scagnostics]{scagnostics}}
#' function from the \code{scagnostics} package.
#' @export
#'
#' @references Wilkinson, Leland, Anushka Anand, and Robert Grossman.
#' "Graph-theoretic scagnostics."
#' Information Visualization, IEEE Symposium on. IEEE Computer Society, 2005
#' @examples
#'  pair_scagnostics(iris)

pair_scagnostics <- function(d, scagnostic = c("Outlying","Skewed","Clumpy","Sparse","Striated",
                                        "Convex","Skinny","Stringy","Monotonic"),
                      handle.na = TRUE, ...) {
  if (!requireNamespace("scagnostics", quietly = TRUE))
    stop("Please install package 'scagnostics' to use pair_scagnostics", call.=FALSE)
  check_df(d)
  scag_choices <- c("Outlying","Skewed","Clumpy","Sparse","Striated",
                    "Convex","Skinny","Stringy","Monotonic")
  sel_scag <- match.arg(scagnostic, scag_choices, several.ok=TRUE)
  
  
  d <- dplyr::select(d, dplyr::where(is.numeric))
  scag1 <- pairwise(d, score = sel_scag[1], pair_type = "nn")
  scag <- scag1
  for (s in sel_scag[-1])
    scag <- rbind(scag, pairwise(d, score = s, pair_type = "nn"))
  scag_fn <- function(x,y) {
    
    x <- d[[x]]
    y <- d[[y]]
    if(handle.na){
      pick <- complete.cases(x, y)
      x <- x[pick]
      y <- y[pick]
    }
    
    scags <- scagnostics::scagnostics(x,y)
    scags[sel_scag]
  }
  
  scag$value <- as.numeric(t(mapply(scag_fn, scag1$x,scag1$y)))
  scag |> 
    dplyr::arrange(.data$x, .data$y)
}
