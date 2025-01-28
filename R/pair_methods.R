
#' Pearson, Spearman or Kendall correlation
#'
#' Calculates one of either pearson, spearman or kendall correlation for every numeric variable pair in a dataset.
#'
#' @param d A dataframe
#' @param method A character string for the correlation coefficient to be calculated. Either "pearson" (default),
#'               "spearman", or "kendall". If the value is "all", then all three correlations are calculated.
#' @param handle.na If TRUE uses pairwise complete observations to calculate correlation coefficient, otherwise NAs not handled.
#' @param ... other arguments
#'
#' @return A tibble of class `pairwise` with calculated association value for every numeric variable pair, 
#' or NULL if there are not at least two numeric variables
#' @seealso See \code{\link{pair_methods}}  for other score options.
#' @export
#'
#' @examples
#' pair_cor(iris)
#' pair_cor(iris, method="kendall")
#' pair_cor(iris, method="spearman")
#' pair_cor(iris, method="all")


pair_cor <- function(d, method="pearson", handle.na=TRUE,...){
  check_df(d)
  d <- d[, sapply(d, is.numeric), drop=FALSE]
  if (ncol(d) > 1){
    if (method == "all"){
      pearson <- pair_cor(d,method="pearson", handle.na = handle.na)
      spearman <- pair_cor(d,method="spearman", handle.na = handle.na)
      kendall <- pair_cor(d,method="kendall", handle.na = handle.na)
      dplyr::bind_rows(pearson ,spearman, kendall)
    }
    else {if (handle.na)
      dcor <- cor(d,method=method,use="pairwise.complete.obs")
    else dcor <- cor(d,method=method,...)
    pairwise(dcor, score=method, pair_type = "nn")
    }
  }
}




ccor <- function(x,y, handle.na=TRUE){
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
  
  d <- d[, sapply(d, function(x) is.numeric(x) | is.factor(x)), drop=FALSE]
  
  if (ncol(d) > 1){
    a <- pairwise(d, score="cancor")
    fn <- function(x,y){
      x <- d[[x]]
      y <- d[[y]]
      ccor(x,y,handle.na=handle.na)
    }
    
    a$value <- mapply(fn, a$x,a$y, USE.NAMES = FALSE)
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
  
  d <- d[, sapply(d, is.numeric), drop=FALSE]
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
    dcor$value <-  mapply(fn, dcor$x,dcor$y, USE.NAMES = FALSE)
    dcor
  }

}




# pair_mine <- function(d, method="mic",handle.na=TRUE,...){
#   if (!requireNamespace("minerva", quietly = TRUE))
#     stop("Please install package 'minerva' to use pair_mine", call.=FALSE)
#   check_df(d)
#   
#   d <- d[, sapply(d, is.numeric), drop=FALSE]
#   if(ncol(d)>1){
#     if (handle.na)
#       dcor <- minerva::mine(d,use="pairwise.complete.obs",...)
#     else dcor <- minerva::mine(d,...)
# 
#     dcor <- dcor[[toupper(method)]]
#     pairwise(dcor, score=method, pair_type = "nn")
#   }
# 
# }


#' MINE family values
#'
#' Calculates MINE family values for every numeric variable pair in a dataset.
#'
#' @param d A dataframe
#' @param method character vector for the MINE value to be calculated. Subset of "MIC","MAS","MEV","MCN","MICR2", "GMIC",  "TIC"  
#' @param handle.na If TRUE uses pairwise complete observations to calculate score, otherwise NAs not handled.
#' @param ... other arguments
#'
#' @return A tibble of class `pairwise` with scores for numeric variable pairs,
#' or NULL if there are not at least two numeric variables
#' @export
#' @details The values are calculated using \code{\link[minerva]{mine}} from \code{minerva}
#' @examples
#'  pair_mine(iris)
#'  pair_mine(iris, method="MAS")

#' @references Reshef, David N., et al. "Detecting novel associations in large data sets."
#' science 334.6062 (2011): 1518-1524
#' 
pair_mine <- function(d, method="MIC",handle.na=TRUE,...){
  if (!requireNamespace("minerva", quietly = TRUE))
    stop("Please install package 'minerva' to use pair_mine", call.=FALSE)
  check_df(d)
  
  d <- d[, sapply(d, is.numeric), drop=FALSE]
  if(ncol(d)>1){
    if (handle.na)
      dcor <- minerva::mine(d,use="pairwise.complete.obs",normalization=TRUE,...)
    else dcor <- minerva::mine(d,normalization=TRUE,...)
    mine_choices <- names(dcor)
    sel_mine <- match.arg(toupper(method), mine_choices, several.ok=TRUE)
    p <- pairwise(dcor[[sel_mine[1]]], score=sel_mine[1], pair_type = "nn")
    for (m in sel_mine[-1]){
      pm <- pairwise(dcor[[m]], score=m, pair_type = "nn")
      p <- rbind(p, pm)
    }
    p |> 
      dplyr::arrange(.data$x, .data$y)
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
#' if (requireNamespace("linkspotter", quietly = TRUE)) { 
#'    pair_nmi(iris)
#' }

pair_nmi <- function(d,handle.na=TRUE,...){
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
      mi <- suppressMessages(linkspotter::maxNMI(x,y))
    }
     nmi$value <-  mapply(fn, nmi$x,nmi$y, USE.NAMES = FALSE)
    nmi
    
  }
}


#' Polychoric correlation
#'
#' Calculates Polychoric correlation using  for every factor variable pair in a dataset.
#'
#' @param d A dataframe
#' @param handle.na ignored. Pairwise complete observations are used automatically.
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
#' pair_polychor(iris)


pair_polychor <- function(d,handle.na=TRUE,...){
  # if (!requireNamespace("polychor", quietly = TRUE))
  #   stop("Please install package 'polychor' to use pair_polychor", call.=FALSE)
  check_df(d)
  
  # polychor automatically does pairwise omit
  d <- d[, sapply(d, is.factor), drop=FALSE]
  if(ncol(d)>1){
    pcor <- pairwise(d, score="polychor", pair_type = "ff")
    pcor$value <- mapply(function(x,y) polycor::polychor(d[[x]],d[[y]],...), pcor$x,pcor$y, USE.NAMES = FALSE)
    pcor
  }

}


#' Polyserial correlation
#'
#' Calculates Polyserial correlation using  for every factor-numeric variable pair in a dataset.
#'
#' @param d A dataframe
#' @param handle.na ignored. Pairwise complete observations are used automatically.
#' @param ... other arguments
#'
#' @return A tibble of class `pairwise` with polyserial correlation for factor-numeric pairs, 
#' or NULL if there are not at least one such pair.
#'
#' @details The polyserial correlation is calculated using the \code{\link[polycor]{polyserial}} function from the
#' \code{polycor} package, and assumes factor levels are in the given order. NAs are automatically handled by pairwise omit.
#' @export
#'
#' @examples
#' pair_polyserial(iris)


pair_polyserial <- function(d,handle.na=TRUE,...){
  # if (!requireNamespace("polychor", quietly = TRUE))
  #   stop("Please install package 'polychor' to use pair_polyserial", call.=FALSE)
  check_df(d)
  
  # polychor automatically does pairwise omit
  d <- dplyr::select(d, dplyr::where(is.numeric) | dplyr::where(is.factor) )
  if(ncol(d)>1){
    pser <- pairwise(d, score="polyserial")
    
    fn <- function(x,y){
      if (is.factor(d[[x]]))
        polycor::polyserial(d[[y]],d[[x]])
      else polycor::polyserial(d[[x]],d[[y]])
    }
     pser <- pser[pser$pair_type == "fn",]
    # if (nrow(pser) == 1) return(NULL)
    suppressWarnings(pser$value <- mapply(fn, pser$x,pser$y, USE.NAMES = FALSE))
    pser
  }
  
}






pair_tau <- function(d,method="B",handle.na=TRUE,...){
  # automatically does pairwise omit for A,B,C
  d <- dplyr::select(d, dplyr::where(is.factor))
  if(ncol(d)>1){
      a <- pairwise(d, score=paste0("tau", method), pair_type = "ff")
      fns <- c("A"= DescTools::KendallTauA, "B"=DescTools::KendallTauB, "C" = DescTools::StuartTauC, "W"=
                 DescTools::KendallW)
      fn <- fns[[method]]
      fnlocal <- function(x,y){
        if (length(unique(d[[x]])) <= 1) return(NA)
        if (length(unique(d[[y]])) <= 1) return(NA)
        if (method =="W")
          # fn(d[c(x,y)], correct=TRUE, na.rm=handle.na,...)
          # version DescTools 0.99.59 compatibility
          fn(cbind(as.integer(d[[x]]),as.integer(d[[y]])), correct=TRUE,...)
        else fn(d[[x]],d[[y]],...)
      }
      a$value <- mapply(fnlocal, a$x,a$y, USE.NAMES = FALSE)
      a
    }
}

#' Kendall's tau  B for association between ordinal factors.
#'
#' Calculates Kendall's tau B every  factor variable pair in a dataset.
#'
#' @param d A dataframe
#' @param ... other arguments
#' @param handle.na ignored. Pairwise complete observations are used automatically.
#' @return  A tibble of class `pairwise` with factor pairs, or NULL if there are not at least two factor variables
#'
#' @details Calculated using  \code{\link[DescTools]{KendallTauB}}. Assumes factor levels are in the given order. 
#' NAs are automatically handled by pairwise omit.
#' @export
#' @examples

#'  d <- data.frame(x=rnorm(20), 
#'                  y=factor(sample(3,20, replace=TRUE)), 
#'                  z=factor(sample(2,20, replace=TRUE)))
#'  pair_tauB(d)

pair_tauB <- function(d,handle.na=TRUE,...){
  if (!requireNamespace("DescTools", quietly = TRUE))
    stop("Please install package 'DescTools' to use pair_tauB ", call.=FALSE)
  check_df(d)
  pair_tau(d, method="B",...)
}

#' Kendall's tau  A for association between ordinal factors.
#'
#' Calculates Kendall's tau A for every  factor variable pair in a dataset.
#'
#' @param d A dataframe
#' @param handle.na ignored. Pairwise complete observations are used automatically.
#' @param ... other arguments
#'
#' @return  A tibble of class `pairwise` with factor pairs, or NULL if there are not at least two factor variables
#'
#' @details Calculated using  \code{\link[DescTools]{KendallTauA}}. Assumes factor levels are in the given order. 
#' NAs are automatically handled by pairwise omit.
#' @export
#' @examples
#'  d <- data.frame(x=rnorm(20), 
#'                  y=factor(sample(3,20, replace=TRUE)), 
#'                  z=factor(sample(2,20, replace=TRUE)))
#'  pair_tauA(d)

pair_tauA <- function(d,handle.na=TRUE,...){
  if (!requireNamespace("DescTools", quietly = TRUE))
    stop("Please install package 'DescTools' to use pair_tauA ", call.=FALSE)
  check_df(d)
  pair_tau(d, method="A",...)
}

#' Stuarts's tau  C for association between ordinal factors.
#'
#' Calculates Stuarts's tau C every  factor variable pair in a dataset.
#'
#' @param d A dataframe
#' @param handle.na ignored. Pairwise complete observations are used automatically.
#' @param ... other arguments
#'
#' @return  A tibble of class `pairwise` with factor pairs, or NULL if there are not at least two factor variables
#'
#' @details Calculated using  \code{\link[DescTools]{StuartTauC}}. Assumes factor levels are in the given order. 
#' NAs are automatically handled by pairwise omit.
#' @export
#' @examples
#'  d <- data.frame(x=rnorm(20), 
#'                  y=factor(sample(3,20, replace=TRUE)), 
#'                  z=factor(sample(2,20, replace=TRUE)))
#'  pair_tauC(d)

pair_tauC <- function(d,handle.na=TRUE,...){
  if (!requireNamespace("DescTools", quietly = TRUE))
    stop("Please install package 'DescTools' to use pair_tauC ", call.=FALSE)
  check_df(d)
  pair_tau(d, method="C")
}

#' Kendall's W for association between ordinal factors.
#'
#' Calculates Kendall's tau W every  factor variable pair in a dataset.
#'
#' @param d A dataframe
#' @param handle.na ignored. Pairwise complete observations are used automatically.
#' @param ... other arguments
#'
#' @return  A tibble of class `pairwise` with factor pairs, or NULL if there are not at least two factor variables
#'
#' @details Calculated using  \code{\link[DescTools]{KendallW}}. Assumes factor levels are in the given order.
#'  NAs are automatically handled by pairwise omit.
#' @export
#' @examples
#'  d <- data.frame(x=rnorm(20), 
#'                  y=factor(sample(3,20, replace=TRUE)), 
#'                  z=factor(sample(2,20, replace=TRUE)))
#'  pair_tauW(d)

pair_tauW <- function(d,handle.na=TRUE,...){
  if (!requireNamespace("DescTools", quietly = TRUE))
    stop("Please install package 'DescTools' to use pair_tauW ", call.=FALSE)
  check_df(d)
  pair_tau(d, method="W")
}

#' Uncertainty coefficient for association between factors.
#'
#' Calculates uncertainty coefficient for every factor variable pair in a dataset.
#'
#' @param d A dataframe
#' @param handle.na ignored. Pairwise complete observations are used automatically.
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
    a$value <- mapply(function(x,y) DescTools::UncertCoef(d[[x]],d[[y]],...), a$x,a$y, USE.NAMES = FALSE)
    a
  }

}


#' Goodman Kruskal's Tau for association between ordinal factors.
#'
#' Calculates Goodman Kruskal's Tau coefficient for every factor variable pair in a dataset.
#'
#' @param d A dataframe
#' @param handle.na ignored. Pairwise complete observations are used automatically.
#' @param ... other arguments
#'
#' @return A tibble of class `pairwise` with Goodman Kruskal's Tau for every factor variable pair,
#' or NULL if there are not at least two factor variables
#' @details The Goodman Kruskal's Tau coefficient is calculated using \code{\link[DescTools]{GoodmanKruskalTau}}
#' function from the \code{DescTools} package. Assumes factor levels are in the given order.
#' NAs are automatically handled by pairwise omit.
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
    fnlocal <- function(x,y) 
      max(DescTools::GoodmanKruskalTau(d[[x]],d[[y]],...),DescTools::GoodmanKruskalTau(d[[y]],d[[x]],...))
    a$value <- mapply(fnlocal, a$x,a$y, USE.NAMES = FALSE)
    a
  }

}


#' Goodman Kruskal's Gamma for association between ordinal factors.
#'
#' Calculates Goodman Kruskal's Gamma coefficient for every factor variable pair in a dataset.
#'
#' @param d A dataframe
#' @param handle.na ignored. Pairwise complete observations are used automatically.
#' @param ... other arguments
#'
#' @return A tibble of class `pairwise` with factor variable pairs and Goodman Kruskal's Gamma coefficient, 
#' or NULL if there are not at least two factor variables

#' @details The Goodman Kruskal's Gamma coefficient is calculated using \code{\link[DescTools]{GoodmanKruskalGamma}}
#' function from the \code{DescTools} package. Assumes factor levels are in the given order.
#'  NAs are automatically handled by pairwise omit.
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
    a$value <- mapply(function(x,y) DescTools::GoodmanKruskalGamma(d[[x]],d[[y]],...), a$x,a$y, USE.NAMES = FALSE)
    a
  }

}

#' Pearson's Contingency Coefficient for association between  factors.
#'
#' Calculates Pearson's Contingency coefficient for every factor variable pair in a dataset.
#'
#' @param d A dataframe
#' @param handle.na ignored. Pairwise complete observations are used automatically.
#' @param ... other arguments
#'
#' @return A tibble of class `pairwise` with calculated Pearson's contingency coefficient for every factor variable
#' pair, or NULL if there are not at least two factor variables
#' @export
#' @details The Pearson's contingency coefficient is calculated using \code{\link[DescTools]{ContCoef}}.
#' NAs are automatically handled by pairwise omit.
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
    a$value <- mapply(function(x,y) DescTools::ContCoef(d[[x]],d[[y]],...), a$x,a$y, USE.NAMES = FALSE)
    a
  }

}


#' Calculates ace based transformations and correlation, handling missing values and factors.
#'
#' @param x a numeric vector or factor
#' @param y a numeric vector or factor
#' @param handle.na  If TRUE uses pairwise complete observations.
#'
#' @return result of acepack::ace
#' @export
#'
#' @examples ace_cor(iris$Sepal.Length, iris$Species)
ace_cor <- function(x,y,handle.na=TRUE) {
  if (!requireNamespace("acepack", quietly = TRUE))
    stop("Please install package 'acepack' to use ace_cor", call.=FALSE)
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
  acepack::ace(x,y, cat=cat)
  
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
#' pair_ace(iris)
pair_ace <- function(d, handle.na = TRUE, ...) {
  if (!requireNamespace("acepack", quietly = TRUE))
    stop("Please install package 'acepack' to use pair_acepack", call.=FALSE)
  check_df(d)
  d <- dplyr::select(d, dplyr::where(is.numeric) | dplyr::where(is.factor) )
  if (ncol(d) > 1){
    ace_score <- pairwise(d, score = "ace")
    ace_fn <- function(x,y) {
      x <- d[[x]]
      y <- d[[y]]
      s <- max(0,ace_cor(x,y,handle.na=handle.na)[["rsq"]]) # sometimes acepack gives very small negatives
      sqrt(s)
    }

    ace_score$value <- mapply(ace_fn, ace_score$x,ace_score$y, USE.NAMES = FALSE)
    ace_score
  }
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
#'  pair_methods

pair_methods <- dplyr::tribble(
  ~name, ~nn, ~ff,  ~fn, ~from, ~range, ~ordinal,
  "pair_cor", TRUE, FALSE, FALSE, "cor", "[-1,1]", NA,
  "pair_dcor", TRUE, FALSE,  FALSE, "energy::dcor2d", "[0,1]", NA,
  "pair_mine", TRUE, FALSE, FALSE, "minerva::mine", "[0,1]",NA,
  "pair_ace", TRUE, TRUE,  TRUE, "acepack::ace", "[0,1]",FALSE,
  "pair_cancor", TRUE, TRUE, TRUE, "cancor", "[0,1]",FALSE,
  "pair_nmi",  TRUE, TRUE,  TRUE, "linkspotter::maxNMI", "[0,1]",FALSE,
  "pair_polychor", FALSE, TRUE, FALSE, "polycor::polychor", "[-1,1]",TRUE,
  "pair_polyserial", FALSE, FALSE, TRUE, "polycor::polyserial", "[-1,1]",TRUE,
  "pair_tauB", FALSE, TRUE, FALSE, "DescTools::KendalTauB", "[-1,1]",TRUE,
  "pair_tauA", FALSE, TRUE, FALSE, "DescTools::KendalTauA", "[-1,1]",TRUE,
  "pair_tauC", FALSE, TRUE, FALSE, "DescTools::StuartTauC", "[-1,1]",TRUE,
  "pair_tauW", FALSE, TRUE, FALSE, "DescTools::KendalW", "[-1,1]",TRUE,
  "pair_gkGamma", FALSE,  TRUE, FALSE, "DescTools::GoodmanKruskalGamma", "[-1,1]",TRUE,
  "pair_gkTau", FALSE,  TRUE, FALSE, "DescTools::GoodmanKruskalTau", "[0,1]",TRUE,
  "pair_uncertainty", FALSE,  TRUE, FALSE, "DescTools::UncertCoef", "[0,1]",FALSE,
  "pair_chi", FALSE,  TRUE,  FALSE, "DescTools::ContCoef", "[0,1]",FALSE,
  "pair_scag", TRUE, FALSE,  FALSE, "scagnostics::scagnostics", "[0,1]",NA,
)


#' Constructs a pairwise result for each level of a by variable.
#'
#' @param d a dataframe
#' @param by a character string for the name of the conditioning variable.
#' @param pair_fun A function returning a `pairwise` from a dataset.
#' @param ungrouped If TRUE calculates the ungrouped score in addition to grouped scores.
#'
#' @return tibble of class "pairwise"
#' @export
#'
#' @examples
#' pairwise_by(iris, by="Species", pair_cor)
#' 
#' 

pairwise_by <- function(d, by, pair_fun, ungrouped=TRUE){
  if (!(by %in% names(d))) cli::cli_abort(c("{.var by} not present in dataset."))
  tab <- table(d[[by]])
  if (any(tab == 1)) cli::cli_abort(c("{by} cannot be used as a grouping variable. Need more than one observation at each level."))
  result <- d |>
    dplyr::group_by(.data[[by]]) |>
    dplyr::group_modify(function(x,y) pair_fun(x)) |>
    dplyr::ungroup() |>
    dplyr::mutate(group=.data[[by]]) |>
    dplyr::select(-dplyr::all_of(by))
  class(result)<-append("pairwise", class(result))
  if (ungrouped){
    overall <- d |>
      dplyr::select(-dplyr::all_of(by)) |>
      pair_fun()
    result <- rbind(result, overall)
  }
  result |> dplyr::arrange(.data$x, .data$y)
}