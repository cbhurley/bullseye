#' A generic function to create a data structure for every variable pair in a dataset
#'
#' Creates a data structure for every variable pair in a dataset.
#'
#' @param data A dataframe.
#' @param score a character string indicating the value of association.
#' @param pair_type a character string specifying the type of variable pair.
#' @return A tbl_df of class `pairwise` for pairs of variables with a column `value` for the score value,
#' `score` for a type of association value and `pair_type` for the type of variable pair.
#' @export
#'
#' @examples
#' pairwise(cor(iris[,1:4]), score="pearson")
#' pairwise(iris)


pairwise <- function(data, score="?", pair_type="?"){
  UseMethod("pairwise", data)
}


#' @describeIn pairwise  pairwise method
#' @export
pairwise.matrix <- function(data, score="?", pair_type="?"){
  m <- data
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
pairwise.data.frame <- function(data, score=NA_character_, pair_type=NA_character_){
  d <- data
  dcor <- diag(ncol(d)) 
  dcor[]<- NA
  rownames(dcor)<- colnames(dcor) <- names(d)
  dcor <- pairwise(dcor, score=score, pair_type=pair_type)
  dcor
}



#' Pearson, Spearman or Kendall correlation
#'
#' Calculates one of either pearson, spearman or kendall correlation for every numeric variable pair in a dataset.
#'
#' @param d A dataframe
#' @param method A character string for the correlation coefficient to be calculated. Either "pearson" (default),
#'               "spearman", or "kendall"
#' @param handle.na If TRUE uses pairwise complete observations to calculate correlation coefficient, otherwise NAs not handled.
#' @param ... other arguments
#'
#' @return A tibble of class `pairwise` with calculated association value for every numeric variable pair, 
#' or NULL if there are not at least two numeric variables
#' @seealso See \code{\link{pairwise_methods}}  for other score options.
#' @export
#'
#' @examples
#' pair_cor(iris)
#' pair_cor(iris, method="kendall")
#' pair_cor(iris, method="spearman")


pair_cor <- function(d, method="pearson", handle.na=TRUE,...){
  check_df(d)
  d <- dplyr::select(d, dplyr::where(is.numeric))
  if (ncol(d) > 1){
    if (handle.na)
      dcor <- cor(d,method=method,use="pairwise.complete.obs")
    else dcor <- cor(d,method=method,...)
    pairwise(dcor, score=method, pair_type = "nn")
  }
}


#' Canonical correlation
#'
#' Calculates canonical correlation for every variable pair in a dataset.
#'
#' @param d A dataframe
#' @param handle.na If TRUE uses pairwise complete observations to calculate correlation coefficient,, otherwise NAs not handled.
#' @param ... other arguments
#'
#' @return A tibble of class `pairwise` with canonical correlation for every numeric or factor or mixed variable pair
#' @export
#'
#' @examples
#' pair_cancor(iris)

pair_cancor <- function(d,handle.na=TRUE,...){
  check_df(d)
  # d <- dplyr::select(d, dplyr::where(is.numeric) | dplyr::where(is.factor) |  dplyr::where(is.character)) |>
  #   dplyr::mutate(dplyr::across(dplyr::where(is.character), as.factor))
  d <- dplyr::select(d, dplyr::where(is.numeric) | dplyr::where(is.factor) )
  if (ncol(d) > 1){
    a <- pairwise(d, score="cancor")
    fn <- function(x,y){
      if(handle.na){
        pick <- complete.cases(x, y)
        x <- x[pick]
        y <- y[pick]
      }
      
      if (length(x) <= 2) {
        return (NA)
      }
      if (!is.numeric(x))
        x <- sapply(unique(x), function(u) as.numeric(x ==u))[,-1]
      if (!is.numeric(y))
        y <- sapply(unique(y), function(u) as.numeric(y ==u))[,-1]
      tryCatch(cancor(x,y)$cor[1], error = function(e) {
        # message("Cannot calculate cancor, returning NA")
        NA
      }
      )
    }
    fn_pair_type <- function(x,y){
      if(is.numeric(x) & is.numeric(y)) {
        "nn"
      } else if(is.factor(x) & is.factor(y)){
        "ff"
      } else {
        "fn"
      }
    }
    a$value <- mapply(function(x,y) fn(d[[x]],d[[y]]), a$x,a$y)
    a$pair_type <- mapply(function(x,y) fn_pair_type(d[[x]],d[[y]]), a$x,a$y)
    a
  }
}

#' Distance correlation
#'
#' Calculates distance correlation for every numeric variable pair in a dataset.
#'
#' @param d A dataframe
#' @param handle.na If TRUE uses pairwise complete observations to calculate distance correlation, otherwise NAs not handled.
#' @param ... other arguments
#' @details The distance correlation is calculated using \code{\link[energy]{dcor2d}} from \code{energy} package
#' @return A tibble of class `pairwise` with distance correlation for every numeric variable pair,
#' or NULL if there are not at least two numeric variables
#' @export
#'
#' @examples
#' pair_dcor(iris)

pair_dcor <- function(d, handle.na=TRUE,...){
  if (!requireNamespace("energy", quietly = TRUE))
    stop("Please install package 'energy' to use pair_dcor", call.=FALSE)
  check_df(d)
  
  d <- dplyr::select(d, dplyr::where(is.numeric))
  if(ncol(d)>1){
    dcor <- pairwise(d, score="dcor", pair_type = "nn")
    fn <- function(x,y){
      x <- d[[x]]
      y <- d[[y]]
      if(handle.na){
        pick <- complete.cases(x, y)
        x <- x[pick]
        y <- y[pick]
      }
      sqrt(energy::dcor2d(x,y,...))
    }
    dcor$value <-  mapply(fn, dcor$x,dcor$y)
    dcor
  }

}

# Easystats correlations
#
# Calculates one of the many correlation coefficients available with easystats package
# for variable pairs in a dataset.
#
# @param d A dataframe
# @param method A character string for the correlation coefficient to be calculated. One of "pearson" (default),
#               "spearman", "kendall", "biserial", "polychoric", "tetrachoric", "biweight", "distance",
#               "percentage" (for percentage bend correlation), "blomqvist" (for Blomqvist's coefficient),
#               "hoeffding" (for Hoeffding's D), "gamma", "gaussian" (for Gaussian Rank correlation) or
#               "shepherd" (for Shepherd's Pi correlation). Setting "auto" will select the most most relevant
#               method depending on the variable types in the dataset.
# @param handle.na NA handling not available
# @param ... other arguments
#
# @return A tibble with calculated association values

#
# @examples
# paireasy(iris)
# paireasy(iris,method="hoeffding")

# paireasy <-function(d,method = "pearson", handle.na=TRUE,...){
#   # no NA handling
#   a <- pairwise(d, score=paste0("EZ", method))
#   ez <- correlation::correlation(d, method=method, ...)[,1:3]
#   #ez <- correlation::correlation(d, method=method)[,1:3]
#   class(ez) <- "data.frame"
#   class(a) <- class(a)[-1]
#   names(ez) <- c("y","x","value")
#   a<-dplyr::rows_patch(a,ez,  by = c("x","y"))
#   class(a) <- append("pairwise",class(a))
#   a
#  }

#' MINE family values
#'
#' Calculates MINE family values for every numeric variable pair in a dataset.
#'
#' @param d A dataframe
#' @param method character string for the MINE value to be calculated. Either "mic" (default), "mas", "mev",
#'               "mcn", or "mic-r2"
#' @param handle.na If TRUE uses pairwise complete observations to calculate score, otherwise NAs not handled.
#' @param ... other arguments
#'
#' @return A tibble of class `pairwise` with scores for numeric variable pairs,
#' or NULL if there are not at least two numeric variables
#' @export
#' @details The values are calculated using \code{\link[minerva]{mine}} from \code{minerva}
#' @examples
#'  pair_mine(iris)
#'  pair_mine(iris, method="mas")

#' @references Reshef, David N., et al. "Detecting novel associations in large data sets."
#' science 334.6062 (2011): 1518-1524


pair_mine <- function(d, method="mic",handle.na=TRUE,...){
  if (!requireNamespace("minerva", quietly = TRUE))
    stop("Please install package 'minerva' to use pair_mine", call.=FALSE)
  check_df(d)
  
  d <- dplyr::select(d, dplyr::where(is.numeric))
  if(ncol(d)>1){
    if (handle.na)
      dcor <- minerva::mine(d,use="pairwise.complete.obs",...)
    else dcor <- minerva::mine(d,...)

    dcor <- dcor[[toupper(method)]]
    pairwise(dcor, score=method, pair_type = "nn")
  }

}


#' Normalized mutual information
#'
#' Calculates normalized mutual information for every numeric or factor or mixed variable pair in a dataset.
#'
#' @param d A dataframe
#' @param handle.na If TRUE uses pairwise complete observations to calculate normalized mutual information, otherwise NAs not handled.
#' @param ... other arguments
#'
#' @details The normalized mutual information is calculated using \code{\link[linkspotter]{maxNMI}} from linkpotter package
#' @return A tibble of class `pairwise`
#' @export
#'
#' @examples
#' pair_nmi(iris)


pair_nmi <- function(d,handle.na=T,...){
  if (!requireNamespace("linkspotter", quietly = TRUE))
    stop("Please install package 'linkspotter' to use pair_nmi",call.=FALSE)
  check_df(d)
  d <- dplyr::select(d, dplyr::where(is.numeric) | dplyr::where(is.factor) )
  if (ncol(d) > 1){
    nmi <- pairwise(d, score="nmi")
    fn <- function(x,y){
      x <- d[[x]]
      y <- d[[y]]
      if(handle.na){
        pick <- complete.cases(x, y)
        x <- x[pick]
        y <- y[pick]
      }
      mi <- linkspotter::maxNMI(x,y)
    }
    fn_pair_type <- function(x,y){
      x <- d[[x]]
      y <- d[[y]]
      if(is.numeric(x) & is.numeric(y)) {
        "nn"
      } else if(is.factor(x) & is.factor(y)){
        "ff"
      } else {
        "nf"
      }
    }
    
    nmi$value <-  mapply(fn, nmi$x,nmi$y)
    nmi$pair_type <- mapply(fn_pair_type, nmi$x,nmi$y)
    nmi
    
  }
}


#' Polychoric correlation
#'
#' Calculates Polychoric correlation using  for every factor variable pair in a dataset.
#'
#' @param d A dataframe
#' @param handle.na If TRUE uses pairwise complete observations to calculate polychoric correlation, otherwise NAs not handled.
#' @param ... other arguments
#'
#' @return A tibble of class `pairwise` with polychoric correlation for factor pairs, 
#' or NULL if there are not at least two factor variables
#'
#' @details The polychoric correlation is calculated using the \code{\link[polycor]{polychor}} function from the
#' \code{polycor} package, and assumes factor levels are in the given order. NAs are automatically handled by pairwise omit.
#' @export
#'
#' @examples
#' pair_polycor(iris)


pair_polycor <- function(d,handle.na=TRUE,...){
  if (!requireNamespace("polycor", quietly = TRUE))
    stop("Please install package 'polycor' to use pair_polycor", call.=FALSE)
  check_df(d)
  
  # polycor automatically does pairwise omit
  d <- dplyr::select(d, dplyr::where(is.factor))
  if(ncol(d)>1){
    pcor <- pairwise(d, score="polycor", pair_type = "ff")
    pcor$value <- mapply(function(x,y) polycor::polychor(d[[x]],d[[y]],...), pcor$x,pcor$y)
    pcor
  }

}

#' Kendall's tau A, B, C and Kendall's W
#'
#' Calculates one of either Kendall's tau A, B, C or Kendall's W for every  factor variable pair in a dataset.
#'
#' @param d A dataframe
#' @param method A character string for the correlation coefficient to be calculated. Either "B" (default),
#'               "A", "C" or "W"
#' @param ... other arguments
#'
#' @return  A tibble of class `pairwise` with factor pairs along with one of either Kendall's tau A, B, C or
#' Kendall's W value, or NULL if there are not at least two factor variables
#'
#' @details The association values Kendall's tau A, B, C or Kendall's W are calculated using \code{\link[DescTools]{KendallTauA}},
#' \code{\link[DescTools]{KendallTauB}}, \code{\link[DescTools]{StuartTauC}} or \code{\link[DescTools]{KendallW}} respectively,from the
#' \code{DescTools} package, and assumes factor levels are in the given order. NAs are automatically handled by pairwise omit.
#'
#' @export
#'
#' @examples
#'  pair_tau(iris)
#'  pair_tau(iris, method="A")
#'  pair_tau(iris, method="C")
#'  pair_tau(iris, method="W")



pair_tau <- function(d,method=c("B","A","C","W"),...){
  if (!requireNamespace("DescTools", quietly = TRUE))
    stop("Please install package 'DescTools' to use pair_tau", call.=FALSE)
  check_df(d)
  
  # automatically does pairwise omit, Kendall
  method <- method[1]
  d <- dplyr::select(d, dplyr::where(is.factor))
  if(ncol(d)>1){
    a <- pairwise(d, score=paste0("tau", method), pair_type = "ff")
    fns <- c("A"= DescTools::KendallTauA, "B"=DescTools::KendallTauB, "C" = DescTools::StuartTauC, "W"=
               DescTools::KendallW)
    # fns <- c("A"= DescTools::KendallTauA, "B"=foo, "C" = DescTools::StuartTauC, "W"=
    #            DescTools::KendallW)
    fn <- fns[[method]]
    fnlocal <- function(x,y){
      if (length(unique(d[[x]])) <= 1) return(NA)
      if (length(unique(d[[y]])) <= 1) return(NA)
      if (method =="W")
        fn(d[c(x,y)], correct=TRUE)
      else fn(d[[x]],d[[y]])
    }
    a$value <- mapply(fnlocal, a$x,a$y)
    a
  }

}

#' Uncertainty coefficient
#'
#' Calculates uncertainty coefficient for every factor variable pair in a dataset.
#'
#' @param d A dataframe
#' @param handle.na If TRUE uses pairwise complete observations to calculate uncertainty coefficient, otherwise NAs not handled.
#' @param ... other arguments
#'
#' @return A tibble of class `pairwise` with every factor variable pair and uncertainty coefficient value,
#' or NULL if there are not at least two factor variables

#' @details The Uncertainty coefficient is calculated using \code{\link[DescTools]{UncertCoef}} function from the
#' \code{DescTools} package.
#'
#' @export
#'
#' @examples
#'  pair_uncertainty(iris)

pair_uncertainty <- function(d,handle.na=TRUE,...){
  if (!requireNamespace("DescTools", quietly = TRUE))
    stop("Please install package 'DescTools' to use pair_uncertainty", call.=FALSE)
  check_df(d)
  
  d <- dplyr::select(d, dplyr::where(is.factor))
  if(ncol(d)>1){
    a <- pairwise(d, score="uncertainty", pair_type = "ff")
    a$value <- mapply(function(x,y) DescTools::UncertCoef(d[[x]],d[[y]],...), a$x,a$y)
    a
  }

}


#' Goodman Kruskal's Tau
#'
#' Calculates Goodman Kruskal's Tau coefficient for every factor variable pair in a dataset.
#'
#' @param d A dataframe
#' @param handle.na If TRUE uses pairwise complete observations, otherwise NAs not handled.
#' @param ... other arguments
#'
#' @return A tibble of class `pairwise` with Goodman Kruskal's Tau for every factor variable pair,
#' or NULL if there are not at least two factor variables
#' @details The Goodman Kruskal's Tau coefficient is calculated using \code{\link[DescTools]{GoodmanKruskalTau}}
#' function from the \code{DescTools} package, and assumes factor levels are in the given order.
#' @export
#'
#' @examples
#'  pair_gkTau(iris)

pair_gkTau <- function(d,handle.na=TRUE,...){
  if (!requireNamespace("DescTools", quietly = TRUE))
    stop("Please install package 'DescTools' to use pair_gkTau", call.=FALSE)
  check_df(d)
  
  d <- dplyr::select(d, dplyr::where(is.factor))
  if(ncol(d)>1){
    a <- pairwise(d, score="gkTau", pair_type = "ff")
    fnlocal <- function(x,y) max(DescTools::GoodmanKruskalTau(d[[x]],d[[y]]),DescTools::GoodmanKruskalTau(d[[y]],d[[x]]))
    a$value <- mapply(fnlocal, a$x,a$y)
    a
  }

}


#' Goodman Kruskal's Gamma
#'
#' Calculates Goodman Kruskal's Gamma coefficient for every factor variable pair in a dataset.
#'
#' @param d A dataframe
#' @param handle.na If TRUE uses pairwise complete observations, otherwise NAs not handled.
#' @param ... other arguments
#'
#' @return A tibble of class `pairwise` with factor variable pairs and Goodman Kruskal's Gamma coefficient, 
#' or NULL if there are not at least two factor variables

#' @details The Goodman Kruskal's Gamma coefficient is calculated using \code{\link[DescTools]{GoodmanKruskalGamma}}
#' function from the \code{DescTools} package,and assumes factor levels are in the given order.
#' @export
#'
#' @examples
#'  pair_gkGamma(iris)

pair_gkGamma <- function(d,handle.na=TRUE,...){
  if (!requireNamespace("DescTools", quietly = TRUE))
    stop("Please install package 'DescTools' to use pair_gkGamma", call.=FALSE)
  
  check_df(d)
  
  d <- dplyr::select(d, dplyr::where(is.factor))
  if(ncol(d)>1){
    a <- pairwise(d, score="gkGamma", pair_type = "ff")
    a$value <- mapply(function(x,y) DescTools::GoodmanKruskalGamma(d[[x]],d[[y]],...), a$x,a$y)
    a
  }

}

#' Pearson's Contingency Coefficient
#'
#' Calculates Pearson's Contingency coefficient for every factor variable pair in a dataset.
#'
#' @param d A dataframe
#' @param handle.na If TRUE uses pairwise complete observations.
#' @param ... other arguments
#'
#' @return A tibble of class `pairwise` with calculated Pearson's contingency coefficient for every factor variable
#' pair, or NULL if there are not at least two factor variables
#' @export
#' @details The Pearson's contingency coefficient is calculated using \code{\link[DescTools]{ContCoef}}
#' function from the \code{DescTools} package.
#'
#' @examples
#'  pair_chi(iris)

pair_chi <- function(d,handle.na=TRUE,...){
  if (!requireNamespace("DescTools", quietly = TRUE))
    stop("Please install package 'DescTools' to use pair_chi", call.=FALSE)
  check_df(d)
  
  d <- dplyr::select(d, dplyr::where(is.factor))
  if(ncol(d)>1){
    a <- pairwise(d, score="chi", pair_type = "ff")
    a$value <- mapply(function(x,y) DescTools::ContCoef(d[[x]],d[[y]],...), a$x,a$y)
    a
  }

}




#' Alternating conditional expectations correlation
#'
#' Calculates the maximal correlation coefficient from alternating conditional expectations algorithm for every variable pair in a dataset.
#'
#' @param d A dataframe
#' @param handle.na If TRUE uses pairwise complete observations, otherwise NAs not handled.
#' @param ... other arguments
#'
#' @return A tibble of class `pairwise` with a maximal correlation  from the alternating conditional expectations
#' algorithm for every variable pair
#'
#' @details The maximal correlation is calculated using alternating conditional expectations
#' algorithm which find the transformations of variables such that the squared correlation
#' is maximised. The \code{\link[acepack]{ace}} function from \code{acepack} package is used for the
#' calculation.
#' @export
#' @references Breiman, Leo, and Jerome H. Friedman.
#' "Estimating optimal transformations for multiple regression and correlation."
#' Journal of the American statistical Association 80.391 (1985): 580-598.
#'
#' @examples
#'  pair_ace(iris)

pair_ace <- function(d, handle.na = T, ...) {
  if (!requireNamespace("acepack", quietly = TRUE))
    stop("Please install package 'acepack' to use pair_acepack", call.=FALSE)
  check_df(d)
  d <- dplyr::select(d, dplyr::where(is.numeric) | dplyr::where(is.factor) )
  if (ncol(d) > 1){
    ace_score <- pairwise(d, score = "ace")
    ace_fn <- function(x,y) {
      
      x <- d[[x]]
      y <- d[[y]]
      if(handle.na){
        pick <- complete.cases(x, y)
        x <- x[pick]
        y <- y[pick]
      }
      cat <- NULL
      if (is.factor(x)) {
        x <- as.numeric(x)
        cat <- 1
      }
      if (is.factor(y)) {
        y <- as.numeric(y)
        cat <- c(cat,0)
      }
      sqrt(acepack::ace(x,y, cat=cat)[["rsq"]])
    }
    
    fn_pair_type <- function(x,y){
      x <- d[[x]]
      y <- d[[y]]
      if(is.numeric(x) & is.numeric(y)) {
        "nn"
      } else if(is.factor(x) & is.factor(y)){
        "ff"
      } else {
        "nf"
      }
    }
    
    ace_score$value <- mapply(ace_fn, ace_score$x,ace_score$y)
    ace_score$pair_type <- mapply(fn_pair_type, ace_score$x,ace_score$y)
    ace_score
  }
}


#' Calculates multiple association measures
#'
#' Calculates multiple association measures for every variable pair in a dataset.
#'
#' @param d dataframe
#' @param scores a vector of names of pairwise functions.
#'
#' @param handle.na If TRUE uses pairwise complete observations to calculate pairwise score, otherwise NAs not handled.
#'
#' @return tibble of class "multi_pairwise"
#' @export
#'
#' @examples
#' iris1 <- iris
#' iris1$Sepal.Length <- cut(iris1$Sepal.Length,3)
#' pair_all(iris1)


pair_all <- function(d,scores=c("pair_cor","pair_dcor","pair_mine","pair_ace",
                                    "pair_cancor","pair_nmi", "pair_polycor",    
                                    "pair_tau","pair_gkGamma","pair_gkTau","pair_uncertainty",
                                    "pair_chi"),
                         handle.na=T) {
  check_df(d)
  results <- vector("list", length(scores))
  for (i in 1:length(scores)){
    results[[i]] <- do.call(what = get(scores[i]), args = list(d = d, handle.na = handle.na))
  }
  
  results <- dplyr::bind_rows(results)
  
  if ("pair_cor" %in% scores) {
    spearman <- pair_cor(d,method="spearman", handle.na = handle.na)
    kendall <- pair_cor(d,method="kendall", handle.na = handle.na)
    results <- dplyr::bind_rows(results, spearman, kendall)
  }
  
  if ("pair_tau" %in% scores) {
    taua <- pair_tau(d,method="A", handle.na = handle.na)
    tauc <- pair_tau(d,method="C", handle.na = handle.na)
    tauw <- pair_tau(d,method="W", handle.na = handle.na)
    results <- dplyr::bind_rows(results, taua, tauc, tauw)
  }
  
  results
  
}



#' Pairwise score functions available in the package
#'
#' A tibble of score functions along with the types of variable pairs these functions
#' can be applied to. It also contains information regarding the packages used to calculate scores
#' and the range of the values calculated.
#'
#' @return tibble
#' @export
#' @examples
#'  pairwise_methods

pairwise_methods <- dplyr::tribble(
  ~name, ~nn, ~ff,  ~fn, ~from, ~range,
  "pair_cor", "y", "n", "n", "cor", "[-1,1]",
  "pair_dcor", "y", "n",  "n", "energy::dcor2d", "[0,1]",
  "pair_mine", "y", "n", "n", "minerva::mine", "[0,1]",
  "pair_ace", "y", "y",  "y", "corVis", "[0,1]",
  "pair_cancor", "y", "y", "y", "corVis", "[0,1]",
  "pair_nmi",  "y", "y",  "y", "linkspotter::maxNMI", "[0,1]",
  "pair_polycor", "n", "y", "n", "polycor::polychor", "[-1,1]",
  "pair_tau", "n", "y", "n", "DescTools::KendalTauA,B,C,W", "[-1,1]",
  "pair_gkGamma", "n",  "y", "n", "DescTools::GoodmanKruskalGamma", "[-1,1]",
  "pair_gkTau", "n",  "y", "n", "DescTools::GoodmanKruskalTau", "[0,1]",
  "pair_uncertainty", "n",  "y", "n", "DescTools::UncertCoef", "[0,1]",
  "pair_chi", "n",  "y",  "n", "DescTools::ContCoef", "[0,1]",
  "pair_scag", "y", "n",  "n", "scagnostics::scagnostics", "[0,1]"
)
