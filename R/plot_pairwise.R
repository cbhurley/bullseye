
#' Pairwise plot in a matrix layout
#'
#' Plots multiple pairwise variable scores in a matrix layout.
#'
#' @param scores The scores for the  matrix plot. Either of class `pairwise` or identical in structure to object of class `pairwise`.
#' @param var_order The variable order to be used. The default NULL means variables in are ordered alphabetically. A value of  
#' "seriate_max" means variables are re-ordered to emphasize pairs with maximum abolute scores. A value of "seriate_max_diff" means 
#' variables are re-ordered to emphasize pairs with maximum score differences. Otherwise var_order must be a subset of variables in scores.
#' @param score_limits a numeric vector of length specifying the limits of the scale. 
#' @param inner_width A number between 0 and 1 specifying radius of the inner bullseye.
#' @param center_level Specifies which level of group goes into the innter bullseye. Defaults to "all".
#' @param na.value used for scores with a value of NA
#' @param pal If provided, should name a one of the sequential or diverging palettes from package colorspace. 
#' See [colorspace::hcl_palettes()]. Otherwise defaults to a blue-red scheme.
#' @param interactive defaults to FALSE
#' @return A `girafe` object if interactive==TRUE, otherwise a `ggplot2`.
#' 
#' If  scores has one value for x,y pair, then a filled circle is drawn with fill representing the score value. If there are multiple values for each x,y pair then the filled circle is split into wedges, with the wedge fill representing the values. 
#' If some rows have `group=center_level`, then the glyph is drawn as a bullseye.
#' If scores has a column `n`, then this is used to weight the slices. See the third example below.
#' @examples
#' pair_cor(iris) |> plot_pairwise()
#' pairwise_scores(iris,by="Species") |> plot_pairwise()
#' pairwise_scores(iris[-(1:30),],by="Species", add.nobs=TRUE) |> plot_pairwise()


#' @export

plot_pairwise <- function(scores, var_order="seriate_max", score_limits=NULL, 
                          inner_width=.5,center_level="all",na.value = "grey80",
                          pal ="Blue-Red 3",interactive=FALSE){
  
  check_pairwise(scores)
  prep <- plot_pairwise_prep(scores, score_limits, var_order=var_order)
  scores <- prep$scores
  score_limits <- prep$score_limits
  var_order <- prep$var_order
  score_label <- prep$score_label
  
  if (length(var_order)==1){
    if (grepl("seriate", var_order)) {
      serfn <- if (var_order == "seriate_max_diff") ser_max_diff else ser_max
      m <- pairwise_to_matrix(scores, serfn, default=0)
      o <- suppressMessages(DendSer::dser(stats::as.dist(-m), cost = DendSer::costLPL))
      var_order <- rownames(m)[o]
    }
  }
  scores$x <- factor(scores$x, levels=var_order)
  scores$y <- factor(scores$y, levels=var_order)
  
  scores2 <- scores
  names(scores2)[match("x", names(scores))] <- "y"
  names(scores2)[match("y", names(scores))] <- "x"
  scores <- rbind(scores, scores2)
  
  diag_df <- scores[1:length(var_order),]
  diag_df$x <- diag_df$y <- factor(var_order, levels=var_order)
  diag_df$value <- NA
  diag_df$text <- diag_df$x
  
  if (is.null(center_level) || !(center_level %in% scores$group)) {
    bullseye <- scores
    scores <- scores[FALSE,]
  }
  else {
    bullseye <- scores[scores$group ==center_level,]
    scores <- scores[scores$group !=center_level,]
  }
  
  
  p <- ggplot(diag_df) +
    facet_grid(ggplot2::vars(.data$x), ggplot2::vars(.data$y)) +
    geom_text(data=diag_df,aes(x=0.05,y=.5,label=.data$text),size=3)+
    theme_void()+
    theme(
      panel.background = element_rect(fill="white", color="grey"),
      legend.position = "bottom",
      strip.text.y = element_blank(),
      strip.text.x = element_blank(),
      panel.spacing = unit(0,'pt'),
      aspect.ratio = 1
    )
  
  scoreslocal <- dplyr::mutate(scores, ytemp= .data$n/sum(.data$n),
                               .by=dplyr::all_of(c("x","y")),
                               boundaries = dplyr::case_when(.data$ytemp == 1 ~ NA_character_ , .default="grey50"))
  
  bullseyelocal <- dplyr::mutate(bullseye, ytemp= .data$n/sum(.data$n),.by=dplyr::all_of(c("x","y")),
                                 boundaries = dplyr::case_when(.data$ytemp == 1 ~ NA_character_ , .default="grey50"))
  cr <- coord_radial(theta="y",inner.radius=.5, expand=FALSE)
  cr$default<- TRUE
  p <- p+
    ggiraph::geom_col_interactive(data=scoreslocal,
                                  aes(x=inner_width, y=.data$ytemp, fill=.data$value, color=.data$boundaries,
                                      tooltip=.data$tooltip),
                                  width=1- inner_width,just=0)+
    cr+
    ggiraph::geom_col_interactive(data=bullseyelocal,
                                  aes(x=0,y=.data$ytemp, fill=.data$value, color=.data$boundaries,
                                      tooltip=.data$tooltip), width= inner_width,just=0)
  
  
  p <- p+ 
    {if (pal %in% rownames(colorspace::hcl_palettes("Diverging")) ) 
      colorspace::scale_fill_continuous_diverging(pal,na.value=na.value,limits=score_limits)
    else if (pal %in% rownames(colorspace::hcl_palettes("Sequential")))
      colorspace::scale_fill_continuous_sequential(pal,na.value=na.value,limits=score_limits)
    else scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b",na.value=na.value,limits=score_limits)}+
      scale_color_identity()+ coord_polar(theta="y") +labs(fill = score_label)
  if (interactive) ggiraph::girafe(ggobj=p) else p
  
}

plot_pairwise_prep <- function(scores, score_limits=NULL, var_order=NULL, ignore_n =FALSE, matrix=TRUE){
   
  mscore <- length(unique(scores$score)) >1
  mgroup <- length(unique(scores$group)) >1
  if (mscore)  score_label <- "scores" 
  else {
    score_label <- unique(scores$score)
    if (is.na(score_label)) score_label <-"scores"
  }
  if (mscore & mgroup)
      scores[[score_label]] <- paste(scores$group,scores$score)
   else if (mgroup)
     scores[[score_label]] <- scores$group
   else  scores[[score_label]] <- scores$score
  
  scores$tooltip <- paste(scores[[score_label]],round(scores$value,2), sep="=")
  
  if (!("n" %in% names(scores))) scores$n <- 1
  if (ignore_n) scores$n <- 1
  
  if (is.null(score_limits)) {
    if (all(scores$value >= 0, na.rm=TRUE) & all(scores$value <= 1, na.rm=TRUE))
      score_limits <- c(0,1)
    else if (all(scores$value >= -1, na.rm=TRUE) & all(scores$value <= 1, na.rm=TRUE))
      score_limits <- c(-1,1)
    else {
      score_limits <- range(scores$value, na.rm=TRUE)
      score_limits <- range(labeling::rpretty(score_limits[1], score_limits[2]))
    }
  }
  
  if (isTRUE(matrix)){
  allvars <- unique(c(scores$x, scores$y))
  if (is.null(var_order))
    var_order <- sort(allvars)
  else if ((length(var_order) == 1) && !grepl("seriate", var_order)) {
    if (length(intersect(allvars, var_order) ==0))
      stop("'var_order' must be NULL, 'seriate_max',  'seriate_max_diff' or a subset of the x and y variables in 'scores'")
    else {
      scores <- scores[scores$x %in% var_order & scores$y %in% var_order,]
    }
  }
  list(scores=scores, var_order=var_order, score_limits=score_limits, score_label=score_label)
  }
  else list(scores=scores, score_limits=score_limits,score_label=score_label)
}





#' Pairwise plot in a linear layout
#'
#' Plots the calculated measures of association among different variable pairs for a dataset in a linear layout.
#'
#' @param scores A tibble with the calculated association measures for the  matrix plot.
#' Either of class `pairwise` or identical in structure to object of class `pairwise`.
#' @param pair_order The variable pair order to be used. The default NULL means pairs are in order of their first appearance in `scores`. A value of  
#' "seriate_max" means pairs are in order of  maximum absolute scores. A value of "seriate_max_diff" means 
#' pairs are in order of maximum scores difference.
#' @param geom The geom to be used. Should be "point" or "tile".
#' @param add_lines When geom= "point" is used, should the points be connected by lines? Defaults to FALSE.
#' @param score_limits a numeric vector of length specifying the limits of the scale. 
#' @param na.value used for geom_tile with a value of NA
#' @param pal For geom="title" only. If provided, should name a one of the sequential or diverging palettes from package colorspace. 
#' See [colorspace::hcl_palettes()]. Otherwise defaults to a blue-red scheme.
#' @param interactive defaults to FALSE
#' @return A `girafe` object if interactive==TRUE, otherwise a `ggplot2`.
#' 
#' @examples
#' plot_pairwise_linear(pairwise_scores(iris))
#' plot_pairwise_linear(pairwise_scores(iris,by="Species"))
#' plot_pairwise_linear(pairwise_multi(iris), geom="point")
#' @export
#' 



plot_pairwise_linear <- function(scores,
                              pair_order = "seriate_max",
                              geom = c("point","tile"),
                              add_lines=FALSE,
                              score_limits=NULL, 
                              na.value = "grey80",pal ="Blue-Red 3",
                              interactive=FALSE){
  check_pairwise(scores)
  geom <- match.arg(geom)
  
  prep <- plot_pairwise_prep(scores, score_limits, matrix=FALSE)
  scores <- prep$scores
  score_limits <- prep$score_limits
  score_label <- prep$score_label
  
  if (all(sapply(2:nrow(scores), function(i) scores$x[1] %in% scores[i, c(1,2)])))
    scores$xy <- unlist(sapply(1:nrow(scores), function(i) setdiff(scores[i, c(1,2)], scores$x[1])))
   else if (all(sapply(2:nrow(scores), function(i) scores$y[1] %in% scores[i, c(1,2)])))
    scores$xy <- unlist(sapply(1:nrow(scores), function(i) setdiff(scores[i, c(1,2)], scores$y[1])))
  else scores$xy <- paste0(scores$x, sep=":", scores$y)
  
  
  if (grepl("seriate", pair_order)){
     serfn <- if (pair_order == "seriate_max_diff") ser_max_diff else ser_max
    ord <- dplyr::summarise(scores,
                     n = dplyr::n(),
                     measure= if (.data$n > 1) serfn(.data$value)  else .data$value,
                     .by=dplyr::all_of(c("xy"))) |> dplyr::arrange(dplyr::desc(.data$measure)) |> dplyr::pull(.data$xy)
    scores$xy <- factor(scores$xy, levels=ord)
  } 
  
  mscore <- length(unique(scores[[score_label]])) >1
  if (geom == "tile"){
    levs <- arrange_tiles_x(scores, score_label)
    scores[[score_label]] <- factor(scores[[score_label]], levels=levs)
    
    p <- ggplot(scores) +
       ggiraph::geom_tile_interactive(aes(x=.data[[score_label]],y=.data$xy,fill=.data$value,
                                tooltip=.data$tooltip)) +
      {if (pal %in% rownames(colorspace::hcl_palettes("Diverging")) ) 
        colorspace::scale_fill_continuous_diverging(pal,na.value=na.value,limits=score_limits)
        else if (pal %in% rownames(colorspace::hcl_palettes("Sequential")))
          colorspace::scale_fill_continuous_sequential(pal,na.value=na.value,limits=score_limits)
        else scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b",na.value=na.value,limits=score_limits)}+
      scale_x_discrete(position = "top", breaks=if (!mscore) NULL else waiver(), expand=c(0,0)) +
      scale_y_discrete(limits=rev, expand=c(0,0)) +
      labs(fill = score_label)+
       theme(panel.background = element_rect(fill=na.value, color=na.value),
            axis.title.x = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
                     axis.text = element_text(size = 8),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()
            )
    
  } else {
    mgroup <- length(unique(scores$group)) >1
    if (mgroup) score_label <- "group"
    p <-  ggplot(scores) +
       {if (identical(score_limits, c(-1,1))) geom_hline(yintercept = 0, color="grey40")} +
      ggiraph::geom_point_interactive(aes(x=.data$xy,y=.data$value,colour=.data[[score_label]], 
                   tooltip=.data$tooltip),
                   show.legend = mscore) +
      {if (isTRUE(add_lines)) geom_line(aes(x=.data$xy,y=.data$value,colour=.data[[score_label]], group= .data[[score_label]]),
                                      show.legend = FALSE)}+
      ylim(score_limits[1],score_limits[2]) +
      coord_flip() +scale_x_discrete(limits=rev) +
      labs(y = "scores")
  }
  p <- p+ theme(legend.position="bottom", axis.title.y  = element_blank())
  if (interactive) ggiraph::girafe(ggobj=p) else p
}


arrange_tiles_x <- function(scores, score_label){
  unique(scores) |> 
    tidyr::pivot_wider(id_cols = c("x","y"),names_from = dplyr::all_of(score_label), 
                       values_from = dplyr::all_of("value") ) |> 
    dplyr::select(-(1:2)) -> d
  if (ncol(d)==1)
    names(d)
  else{
    dnames <- names(d)
    ov <- grep("all", dnames)
    dnames[ov] <- paste("000",dnames[ov])
    d <- d[,order(dnames)]
    d1 <- t(0+is.na(as.matrix(d)))
    if (length(unique(scores$group)) > 1) {
    word1 <- sub(" .*", "",rownames(d1))
    d1 <- d1+  match(word1, unique(word1))*10
    }
    o <- DendSer::dser(stats::dist(d1, method="manhattan"), ser_weight = rowSums(d1),cost = DendSer::costLS)
    names(d)[o]
  }
}


#' Plot method for class `pairwise`. 
#'
#' @param x An object of class `pairwise`
#' @param type If "matrix", calls `plot_pairwise`, if "linear" calls `plot_pairwise_linear`
#' @param ... further arguments to \code{plot_pairwise} or \code{plot_pairwise_linear}
#'
#' @return a plot
#' @export
#'
#' @examples
#' plot(pairwise_scores(iris))
plot.pairwise<- function(x, type=c("matrix", "linear"), ...){
  if (type[1]=="matrix") 
    plot_pairwise(x,...)
  else if (type[1]=="linear") 
  plot_pairwise_linear(x,...)
}


ser_max_diff <- function(x){
  if (all(is.na(x))) 0 else diff(range(x, na.rm=TRUE)) 
}

ser_max <- function(x){
  if (all(is.na(x))) 0 else max(abs(x), na.rm=TRUE)
}



#' Converts a pairwise to a symmetric matrix. Uses the first entry for each (x,y) pair.
#' @param x An object of class pairwise
#' @return A symmetric matrix
#' @export