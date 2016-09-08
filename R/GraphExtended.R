# Extended Graph 
# Author: Walter Xie
# Accessed on 15 May 2016

validateAttr <- function(attr.df, colour.id=NULL, shape.id=NULL, link.id=NULL,
                         text.id=NULL, text.size.id=NULL) {
  attr.v <- c()
  if (! is.null(colour.id)) {
    if (! colour.id %in% colnames(attr.df))
      stop("Invalid colour.id,", colour.id,  "not exsit in meta data column names !\n")
    attr.v <- c(attr.v, colour.id)
  } 
  if (! is.null(link.id)) { 
    if (! link.id %in% colnames(attr.df) )
      stop("Invalid link.id,", link.id,  "not exsit in meta data column names !\n")
    attr.v <- c(attr.v, link.id)
  } 
  if (! is.null(shape.id)) { 
    if (! shape.id %in% colnames(attr.df) )
      stop("Invalid shape.id,", shape.id,  "not exsit in meta data column names !\n")
    attr.v <- c(attr.v, shape.id)
  }
  if (! is.null(text.id)) { 
    if (! text.id %in% colnames(attr.df) )
      stop("Invalid shape.id,", text.id,  "not exsit in meta data column names !\n")
    attr.v <- c(attr.v, text.id)
  }
  if (! is.null(text.size.id)) { 
    if (! text.size.id %in% colnames(attr.df) )
      stop("Invalid text.size.id,", text.size.id,  "not exsit in meta data column names !\n")
    attr.v <- c(attr.v, text.size.id)
  }
  return(attr.v)
}

#' @name extScatterPlot
#' @title The Graphs Extended From \code{\link{ggScatterPlot}}
#'
#' @description 
#' The difference between principal components analysis (PCA) 
#' and multidimensional scaling (MDS) is discussed in 
#' \url{http://stats.stackexchange.com/questions/14002/
#' whats-the-difference-between-principal-components-analysis-and-multidimensional}.
#' 
#' A classical MDS is also known as principal coordinates analysis (Gower, 1966).
#' 
#' @note Because functions of computing NMDS \code{\link{metaMDS}}, 
#' PCoA \code{\link{cmdscale}}, and PCA \code{\link{prcomp}} 
#' take different input data, please make sure you give a valid input to plot. 
#' 
#' @return 
#' A \code{\link{gtable}} object, which turns off clipping.
#' It needs to use \code{\link{pdf.gtplot}} to create pdf, 
#' or \code{\link{plot.gtable}} to plot in console. 
#' 
#' @param attr.df A data frame of meta data to define 
#' how the plot is coloured, shaped, or linked.
#' The row names have to match row names of the input data, 
#' such as a matix or \code{\link{dist}} object.
#' The column names have to contain the id of colour, shape, link, or text.
#' 
#' This data frame will \code{\link{merge}} with data points generated 
#' by PCA or MDS, so that the rows whose names are not matched are dropped. 
#' @param colour.id,shape.id,link.id The column name in \code{attr.df} to 
#' define how the data points are coloured, shaped, or linked according to values.
#' @param text.or.point If 1, then display the text only; if 2, then only points; 
#' if 3, the default, then both the text and points.
#' @param k The number of dimensions in \code{\link{metaMDS}}, 
#' or maximum dimension of the space in \code{\link{cmdscale}}. Default to 2.
#' @param ... Other arguments passed to \code{\link{ggScatterPlot}}.

#' @details 
#' NMDS plot \code{gtNMDSPlot} uses \code{\link{metaMDS}} in \code{\link{vegan}}
#' to create Nonmetric Multidimensional Scaling (NMDS) plot.
#' 
#' @param comm Community data for NMDS plot. 
#' It is dissimilarities either as a \code{\link{dist}} structure 
#' or as a community data defined by \code{\link{vegan}} which is 
#' a transposed matrix of community matrix in \code{\link{ComMA}}.
#' Use \code{\link{transposeDF}} to rotate the data frame.
#' Detail to \code{\link{metaMDS}}.
#' @param distance Dissimilarity index used in \code{\link{vegdist}}.
#' Default to "bray".
#' @param title.add.stress If TRUE, the default, then add stress from 
#' \code{\link{metaMDS}} to the title.
#' @keywords graph
#' @export
#' @examples
#' nmds.plot <- gtNMDSPlot(comm, env, colour.id="FishSpecies", shape.id="FeedingPattern", add.text=T)
#' plot(nmds.plot)
#' 
#' @rdname extScatterPlot
gtNMDSPlot <- function(comm, attr.df, colour.id=NULL, shape.id=NULL, link.id=NULL, 
                       text.id=NULL, text.or.point=3, text.size.id=NULL, text.size=3,
                       distance="bray", k = 2, title="NMDS", title.add.stress=TRUE, 
                       verbose=TRUE, ...) {
  if (! missing(attr.df)) {
    if (! all(rownames(as.matrix(comm)) %in% rownames(attr.df)) )
      stop("Invalid attr.df,", paste(rownames(as.matrix(comm)), collapse = ","), 
           "should match", paste(rownames(attr.df), collapse = ","), "!\n")
  }
  
  # Run metaMDS, get points and stress
  suppressMessages(require(vegan))
  mds <- metaMDS(comm, distance = distance, k = k)
  
  if (is.null(mds))
    stop("No points returned from metaMDS, please check your input data !\n")
  
  df.points <- as.data.frame(mds$points)
  df.points <- df.points[order(rownames(df.points)),]
  
  if (title != "" && title.add.stress)
    title <- paste0(title, " (stress ", round(mds$stress, 2), ")")
  
  if (! missing(attr.df)) {
    #rownames(df.points) <- tolower(rownames(df.points))
    #rownames(attr.df) <- tolower(rownames(attr.df))
    
    validateAttr(attr.df, colour.id=colour.id, shape.id=shape.id, link.id=shape.id,
                 text.size.id)
    
    df.points.merge <- merge(df.points, attr.df, by = "row.names")
    
    if (nrow(df.points.merge) != nrow(df.points) || nrow(df.points.merge) != nrow(attr.df)) 
      warning(paste("Some data points are missing after merge ! nrow(df.points.merge) =", 
                    nrow(df.points.merge), ", nrow(df.points) =", nrow(df.points), 
                    ", nrow(attr.df) =", nrow(attr.df) ))
  } else {
    df.points.merge <- df.points
    df.points.merge[,"Row.names"] <- rownames(df.points) 
  }
  
  if (text.or.point != 2)
    text.id="Row.names"
  
  if (! is.null(text.size.id)) {
    # renormalization min to 1
    min.s=min(df.points.merge[,text.size.id])
    max.s=max(df.points.merge[,text.size.id])
    text.size=round(df.points.merge[,text.size.id]/min.s, 2)
  }

  # Plot MDS ordination
  gg.plot <- ComMA::ggScatterPlot(df.points.merge, x.id="MDS1", y.id="MDS2", 
                                  colour.id=colour.id, shape.id=shape.id, link.id=link.id, 
                                  text.id=text.id, text.or.point=text.or.point, 
                                  text.size=text.size, title=title, legend.title.size=text.size.id,
                                  verbose=verbose, ...)
  # turns off clipping
  gt <- ComMA::unclip.ggplot(gg.plot) 
  
  return(gt)
}

#' @details
#' Classical MDS \code{gtClassicalMDSPlot} uses \code{\link{cmdscale}} to  
#' create a classical (metric) multidimensional scaling plot of a data matrix. 
#' Also known as principal coordinates analysis (Gower, 1966).
#' 
#' @param dist.comm A distance structure such as that returned by \code{\link{dist}},  
#' or a full symmetric matrix containing the dissimilarities for the classical 
#' MDS \code{\link{cmdscale}}.
#' @param eig Indicates whether eigenvalues should be returned from \code{\link{cmdscale}}. 
#' Default to TRUE.
#' @param add.const Logical indicating in \code{\link{cmdscale}} if an additive 
#' constant c* should be computed, and added to the non-diagonal dissimilarities 
#' such that the modified dissimilarities are Euclidean. Default to FALSE.
#' @keywords graph
#' @export
#' @examples
#' cmds.plot <- gtClassicalMDSPlot(dist.comm, env, colour.id="FishSpecies", shape.id="FeedingPattern")
#' plot(cmds.plot)
#' 
#' @rdname extScatterPlot
gtClassicalMDSPlot <- function(dist.comm, attr.df, 
                               colour.id=NULL, shape.id=NULL, link.id=NULL, 
                               text.id=NULL, text.or.point=3, 
                               text.size.id=NULL, text.size=3,
                               eig=TRUE, k=2, add.const=FALSE, 
                               title="Classical MDS", verbose=TRUE, ...) {
  if (! missing(attr.df)) {
    if (! all(rownames(as.matrix(dist.comm)) %in% rownames(attr.df)) )
      stop("Invalid attr.df,", paste(rownames(as.matrix(dist.comm)), collapse = ","), 
           "should match", paste(rownames(attr.df), collapse = ","), "!\n")
  }
  
  # Run cmdscale, get points, class(fit) == matrix
  fit <- cmdscale(dist.comm, eig=TRUE, k=2, add = add.const) 
  
  df.points <- as.data.frame(fit$points)
  colnames(df.points) <- c("MDS1", "MDS2")

  if (! missing(attr.df)) {
    #rownames(df.points) <- tolower(rownames(df.points))
    #rownames(attr.df) <- tolower(rownames(attr.df))
    
    validateAttr(attr.df, colour.id=colour.id, shape.id=shape.id, link.id=shape.id)
    
    df.points.merge <- merge(df.points, attr.df, by = "row.names")
    
    if (nrow(df.points.merge) != nrow(df.points) || nrow(df.points.merge) != nrow(attr.df)) 
      warning(paste("Some data points are missing after merge ! nrow(df.points.merge) =", 
                    nrow(df.points.merge), ", nrow(df.points) =", nrow(df.points), 
                    ", nrow(attr.df) =", nrow(attr.df) ))
  } else {
    df.points.merge <- df.points
    df.points.merge[,"Row.names"] <- rownames(df.points) 
  }
  
  if (text.or.point != 2)
    text.id="Row.names"
  
  if (! is.null(text.size.id)) {
    if (! text.size.id %in% colnames(attr.df))
        stop("Invalid text.size.id,", text.size.id,  "not exsit in meta data column names !\n")
    # ggplot2 sizes are measured in mm. 
    # Therefore cex=1 would correspond to size=3.81 or size=5.08 
    # depending on if it is being scaled by the width or height.
    # renormalization min to 1
    min.s=min(df.points.merge[,text.size.id])
    max.s=max(df.points.merge[,text.size.id])
    text.size=round(df.points.merge[,text.size.id]/min.s*3.81, 2)
  }
  
  # Plot MDS ordination
  gg.plot <- ComMA::ggScatterPlot(df.points.merge, x.id="MDS1", y.id="MDS2", 
                                  colour.id=colour.id, shape.id=shape.id, link.id=link.id, 
                                  text.id=text.id, text.or.point=text.or.point, 
                                  text.size=text.size, title=title, 
                                  verbose=verbose, ...)
  # turns off clipping
  gt <- ComMA::unclip.ggplot(gg.plot) 
  
  return(gt)
}



#' @details
#' PCA plot \code{gtPCAPlot} uses \code{\link{prcomp}} to create 
#' a principal components analysis (PCA) plot.
#' 
#' @param df.comm A numeric or complex matrix (or data frame), 
#' which provides the data for the principal components analysis 
#' \code{\link{prcomp}}.
#' For example, a community data defined by \code{\link{vegan}} which is 
#' a transposed matrix of community matrix in \code{\link{ComMA}}.
#' Use \code{\link{transposeDF}} to rotate the data frame.
#' @param scale.pca Logical indicating in \code{\link{prcomp}} 
#' whether the variables should be scaled to have unit variance 
#' before the analysis takes place. Default to TURE.
#' @keywords graph
#' @export
#' @examples
#' pca.plot <- gtPCAPlot(df.comm, env, colour.id="FishSpecies", shape.id="FeedingPattern")
#' plot(pca.plot)
#' 
#' @rdname extScatterPlot
gtPCAPlot <- function(df.comm, attr.df, x.i=1, y.i=2, 
                      colour.id=NULL, shape.id=NULL, link.id=NULL, 
                      text.id=NULL, text.or.point=3, text.data=NULL, 
                      scale.pca=TRUE, title="PCA", 
                      verbose=TRUE, ...) {
  
  if (! missing(attr.df)) {
    if (! all(rownames(as.matrix(comm)) %in% rownames(attr.df)) )
      stop("Invalid attr.df,", paste(rownames(as.matrix(comm)), collapse = ","), 
           "should match", paste(rownames(attr.df), collapse = ","), "!\n")
  }
  
  # Run prcomp, get points 
  pca <- prcomp(comm, scale. = scale.pca)
  df.points <- as.data.frame(pca$rotation)
  df.points <- df.points[order(rownames(df.points)),]
  
  if (x.i>=y.i || x.i < 1 || y.i > ncol(df.points))
    stop("Invalid x.i", x.i,  "or y.i", y.i, "for PCA dimension index !\n")
  
  if (! missing(attr.df)) {
    #rownames(df.points) <- tolower(rownames(df.points))
    #rownames(attr.df) <- tolower(rownames(attr.df))
    
    validateAttr(attr.df, colour.id=colour.id, shape.id=shape.id, link.id=shape.id)
    
    df.points.merge <- merge(df.points, attr.df, by = "row.names")
    
    if (nrow(df.points.merge) != nrow(df.points) || nrow(df.points.merge) != nrow(attr.df)) 
      warning(paste("Some data points are missing after merge !", 
                    nrow(df.points.merge), "!=", nrow(df.points), "!=", nrow(attr.df) ))
  } else {
    df.points.merge <- df.points
    df.points.merge[,"Row.names"] <- rownames(df.points) 
  }
  
  if (text.or.point != 2)
    text.id="Row.names"
  
  x.id <- paste0("PC", x.i)
  y.id <- paste0("PC", y.i)
  # Plot MDS ordination
  gg.plot <- ComMA::ggScatterPlot(df.points.merge, x.id="PC1", y.id="PC2", 
                                  colour.id=colour.id, shape.id=shape.id, link.id=link.id, 
                                  text.id=text.id, text.or.point=text.or.point, 
                                  text.data=text.data, title=title, 
                                  verbose=verbose, ...)
  # turns off clipping
  gt <- ComMA::unclip.ggplot(gg.plot) 
  
  return(gt)
}

#' Rarefaction curves for multi-sample, which is extended from \code{\link{gtLine}}. 
#' 
#' @param df.size A data frame required to \code{\link{merge}} with \code{attr.df} 
#' and \code{\link{melt}} before making a group of rarefaction curves.
#' The rows are samples and must be a subset of rows in \code{attr.df}, 
#' columns are the subsampled data points to draw the curve. 
#' The data frame of rarefaction curve of phylogenetic diversity can be generated by 
#' \code{\link{getPhylorareDF}}. 
#' For example,
#' \tabular{rrrrr}{
#'   Samples \tab size.1 \tab size.5 \tab size.67 \tab ...\cr
#'   Sample1 \tab 1.845445 \tab 3.679956 \tab 9.191672 \tab ...\cr
#'   Sample2 \tab 2.047155 \tab 10.41827 \tab 17.34067 \tab ...\cr
#'   Sample3 \tab 0.06646017 \tab 1.65030905 \tab NaN \tab ... 
#' } 
#' @param attr.df A data frame to define how rarefaction curves are coloured, shaped.
#' @param group.id A column name from \code{df.to.melt} to \code{\link{melt}} 
#' subsampled data points, which is copied from merged "row.names".
#' @param line.or.point If 1, then display the line only; if 2, then only points; 
#' if 3, the default, then both the line and points.
#' @param end.point.only If TRUE, as default, then only display the points in 
#' the end of lines. Otherwise, display all the points.
#' @param x.prefix The regular expression to remove prefix from column names in 
#' \code{df.to.melt} (e.g. size.100). Default to "^.*?\\\\.".
#' @param ... Other arguments passed to \code{\link{gtLine}}.
#' @keywords graph
#' @export
#' @examples 
#' rare.curv <- gtRarefactionCurve(df.size, attr.df, group.id="Samples", colour.id="Species",  
#'                        shape.id="GutSegment",point.size=2, x.trans="log", auto.scale.x=T)
#' plot(rare.curv)
gtRarefactionCurve <- function(df.size, attr.df, group.id="Samples", colour.id=NULL,
                               line.type = 2, line.alpha=0.75, shape.id=NULL,
                               x.prefix="^.*?\\.", line.or.point=3, end.point.only=TRUE,
                               title="Rarefaction Curves", x.lab="Reads", y.lab="Diversity", 
                               verbose=TRUE, ...) {
  if (! missing(attr.df)) {
    if (! all(rownames(df.size) %in% rownames(attr.df)) )
      stop("Invalid attr.df,", paste(rownames(df.size), collapse = ","), 
           "should match", paste(rownames(attr.df), collapse = ","), "!\n")
    
    attr.v <- validateAttr(attr.df, colour.id=colour.id, shape.id=shape.id)
    
    # cannot merge 1-column df
    merge.df <- merge(df.size, attr.df, by = "row.names", sort = FALSE)
    colnames(merge.df)[1] <- group.id
    # only take column needed, make sure not causing problem after melt
    merge.df <- merge.df[,c(group.id, colnames(df.size), attr.v)]
  } else {
    merge.df <- df.size
    merge.df[, group.id] <- rownames(df.size)
  }
  if (verbose)
    cat("colnames(merge.df) = ", paste(colnames(merge.df), collapse = ","), "\n")
  
  suppressMessages(require(reshape2))
  if (verbose)
    cat("melt id = ", paste(c(group.id, attr.v), collapse = ","), "\n")
  # melt by group.id + attr.v
  melt.df <- melt(merge.df, id=c(group.id, attr.v))
  # is.na("NaN") = FALSE
  melt.df[melt.df =="NaN"] <- NA_character_
  # rm all NaN
  melt.df <- melt.df[complete.cases(melt.df),]
  
  # rm all prefix "size."
  if (! is.null(x.prefix))
    melt.df$variable <- gsub(x.prefix, "", melt.df$variable)
  
  melt.df$variable <- as.numeric(melt.df$variable)
  melt.df$value <- as.numeric(melt.df$value)
  
  if (is.na(melt.df$variable) || is.na(melt.df$value))
    stop("melt.df$variable or melt.df$value has NA, please give numberic values !")
  
  point.data=NULL
  if (end.point.only) {
    aggr.formula <- paste("cbind(variable, value) ~ Samples")
    if (! missing(attr.df))
      aggr.formula <- paste(aggr.formula, "+", paste(attr.v, collapse = "+"))
    if (verbose)
      cat("aggregate formula for point.data = (", aggr.formula, ")\n")
    
    point.data <- aggregate(as.formula(aggr.formula), melt.df, max)
  } 
  
  gg.plot <- ComMA::ggLineWithPoints(melt.df, x.id="variable", y.id="value", group.id=group.id, 
                                     colour.id=colour.id, shape.id=shape.id, 
                                     point.data=point.data, line.or.point=line.or.point,
                                     line.type=line.type, line.alpha=line.alpha,
                                     title=title, x.lab=x.lab, y.lab=y.lab, 
                                     verbose=verbose, ...)
  # turns off clipping
  gt <- ComMA::unclip.ggplot(gg.plot) 
  
  return(gt)
}

