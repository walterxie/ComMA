# Extended Graph 
# Author: Walter Xie
# Accessed on 15 May 2016

validateAttr <- function(attr.df, colour.id=NULL, shape.id=NULL, link.id=NULL) {
  attr.v <- c()
  if (! is.null(colour.id)) {
    if (! colour.id %in% colnames(attr.df))
      stop("Invalid colour.id,", colour.id,  "not exsit in column names !\n")
    attr.v <- c(attr.v, colour.id)
  } 
  if (! is.null(link.id)) { 
    if (! link.id %in% colnames(attr.df) )
      stop("Invalid link.id,", link.id,  "not exsit in column names !\n")
    attr.v <- c(attr.v, link.id)
  } 
  if (! is.null(shape.id)) { 
    if (! shape.id %in% colnames(attr.df) )
      stop("Invalid shape.id,", shape.id,  "not exsit in column names !\n")
    attr.v <- c(attr.v, shape.id)
  }
  return(attr.v)
}

#' NMDS Plot
#' 
#' Use \code{\link{metaMDS}} in \code{\link{vegan}}
#' to create Nonmetric Multidimensional Scaling (NMDS) plot.
#' It is extended from \code{\link{gtScatterPlot}}. 
#' 
#' @param comm Community data for NMDS plot. 
#' It is dissimilarities either as a \code{\link{dist}} structure 
#' or as a community data defined by \code{\link{vegan}} which is 
#' a transposed matrix of community matrix in \code{\link{ComMA}}.
#' Use \code{\link{transpose.df}} to rotate the data frame.
#' Detail to \code{\link{metaMDS}}.
#' @param distance Dissimilarity index used in \code{\link{vegdist}}.
#' @param attr.df A data frame to define how NMDS plot is coloured, shaped, or linked.
#' where its row names have to match row names of community data \code{comm}.
#' @param ... Other arguments passed to \code{\link{gtScatterPlot}}.
#' @keywords graph
#' @export
#' @examples
#' 
#' nmds.plot <- gtNMDSPlot(comm, env, colour.id="FishSpecies", shape.id="FeedingPattern", add.text=T)
#' require(grid)
#' grid.draw(nmds.plot)
gtNMDSPlot <- function(comm, attr.df, colour.id=NULL, shape.id=NULL, link.id=NULL, 
                       add.text=TRUE, text.data = NULL, text.size=3, palette=NULL,
                       xintercept=NULL, yintercept=NULL, 
                       distance="bray", title="MDS", verbose=TRUE, ...) {
  if (! missing(attr.df)) {
    if (! all(rownames(as.matrix(comm)) %in% rownames(attr.df)) )
      stop("Invalid attr.df,", paste(rownames(as.matrix(comm)), collapse = ","), 
           "should match", paste(rownames(attr.df), collapse = ","), "!\n")
  }
  
  # Run metaMDS, get points and stress
  suppressMessages(require(vegan))
  mds <- metaMDS(comm, distance = distance)
  pts.mds <- as.data.frame(mds$points)
  pts.mds <- pts.mds[order(rownames(pts.mds)),]
  
  if (title != "")
    title <- paste0(title, " (stress ", round(mds$stress, 2), ")")
  
  if (! missing(attr.df)) {
    #rownames(pts.mds) <- tolower(rownames(pts.mds))
    #rownames(attr.df) <- tolower(rownames(attr.df))
    
    validateAttr(attr.df, colour.id=colour.id, shape.id=shape.id, link.id=shape.id)
    
    pts.mds.merge <- merge(pts.mds, attr.df, by = "row.names")
    
    if (nrow(pts.mds.merge) != nrow(pts.mds) || nrow(pts.mds.merge) != nrow(attr.df)) 
      warning(paste("Some data points are missing after merge ! nrow(pts.mds.merge) =", 
                    nrow(pts.mds.merge), ", nrow(pts.mds) =", nrow(pts.mds), ", nrow(attr.df) =", nrow(attr.df) ))
  } else {
    pts.mds.merge <- pts.mds
    pts.mds.merge[,"Row.names"] <- rownames(pts.mds) 
  }
  
  if (add.text)
    text.id="Row.names"
  
  # Plot MDS ordination
  gt <- ComMA::gtScatterPlot(pts.mds.merge, x.id="MDS1", y.id="MDS2", colour.id=colour.id, 
                             shape.id=shape.id, link.id=link.id, text.id=text.id, 
                             palette=palette, xintercept=xintercept, yintercept=yintercept,
                             text.data=text.data, text.size=text.size, title=title, 
                             verbose=verbose, ...)
  
  return(gt)
}


#' PCA Plot
#' 
#' Use \code{\link{prcomp}} to create Principal Components Analysis (PCA) plot.
#' It is extended from \code{\link{gtScatterPlot}}. 
#' 
#' @param comm Community data for NMDS plot. 
#' It is dissimilarities either as a \code{\link{dist}} structure 
#' or as a community data defined by \code{\link{vegan}} which is 
#' a transposed matrix of community matrix in \code{\link{ComMA}}.
#' Use \code{\link{transpose.df}} to rotate the data frame.
#' Detail to \code{\link{prcomp}}.
#' @param attr.df A data frame to define how NMDS plot is coloured, shaped, or linked.
#' where its row names have to match row names of community data \code{comm}.
#' @param ... Other arguments passed to \code{\link{gtScatterPlot}}.
#' @keywords graph
#' @export
#' @examples
#' 
#' pca.plot <- gtNMDSPlot(comm, env, colour.id="FishSpecies", shape.id="FeedingPattern", add.text=T)
#' require(grid)
#' grid.draw(pca.plot)
gtPCAPlot <- function(comm, attr.df, x.i=1, y.i=2, colour.id=NULL, shape.id=NULL, link.id=NULL, 
                      add.text=TRUE, text.data=NULL, text.size=3, palette=NULL,
                      xintercept=NULL, yintercept=NULL, title="PCA", verbose=TRUE, ...) {
  if (! missing(attr.df)) {
    if (! all(rownames(as.matrix(comm)) %in% rownames(attr.df)) )
      stop("Invalid attr.df,", paste(rownames(as.matrix(comm)), collapse = ","), 
           "should match", paste(rownames(attr.df), collapse = ","), "!\n")
  }
  
  # Run prcomp, get points 
  pca <- prcomp(comm, scale. = TRUE)
  pts.pca <- as.data.frame(pca$rotation)
  pts.pca <- pts.pca[order(rownames(pts.pca)),]
  
  if (x.i>=y.i || x.i < 1 || y.i > ncol(pts.pca))
    stop("Invalid x.i", x.i,  "or y.i", y.i, "for PCA dimension index !\n")
  
  if (! missing(attr.df)) {
    #rownames(pts.pca) <- tolower(rownames(pts.pca))
    #rownames(attr.df) <- tolower(rownames(attr.df))
    
    validateAttr(attr.df, colour.id=colour.id, shape.id=shape.id, link.id=shape.id)
    
    pts.pca.merge <- merge(pts.pca, attr.df, by = "row.names")
    
    if (nrow(pts.pca.merge) != nrow(pts.pca) || nrow(pts.pca.merge) != nrow(attr.df)) 
      warning(paste("Some data points are missing after merge !", 
                    nrow(pts.pca.merge), "!=", nrow(pts.pca), "!=", nrow(attr.df) ))
  } else {
    pts.pca.merge <- pts.pca
    pts.pca.merge[,"Row.names"] <- rownames(pts.pca) 
  }
  
  if (add.text)
    text.id="Row.names"
  
  x.id <- paste0("PC", x.i)
  y.id <- paste0("PC", y.i)
  # Plot MDS ordination
  gt <- ComMA::gtScatterPlot(pts.pca.merge, x.id="PC1", y.id="PC2", colour.id=colour.id, 
                             shape.id=shape.id, link.id=link.id, text.id=text.id, 
                             palette=palette, xintercept=xintercept, yintercept=yintercept,
                             text.data=text.data, text.size=text.size, title=title, 
                             verbose=verbose, ...)
  
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
#' @param x.prefix The regular expression to remove prefix from column names in 
#' \code{df.to.melt} (e.g. size.100). Default to "^.*?\\\\.".
#' @param ... Other arguments passed to \code{\link{gtLine}}.
#' @keywords graph
#' @export
#' @examples 
#' rare.curv <- gtRarefactionCurve(df.size, attr.df, group.id="Samples", colour.id="Species", shape.id="GutSegment", 
#'              point.size=2, x.trans="log", auto.scale.x=T)
#' require(grid)
#' grid.draw(rare.curv)
gtRarefactionCurve <- function(df.size, attr.df, group.id="Samples", colour.id=NULL,
                               x.prefix="^.*?\\.", end.point.only=TRUE,
                               title="Rarefaction Curves", x.lab="Reads", y.lab="Diversity", ...) {
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
  
  gt <- ComMA::gtLine(melt.df, x.id="variable", y.id="value", group.id=group.id, 
                      colour.id=colour.id, shape.id=shape.id, 
                      title=title, x.lab=x.lab, y.lab=y.lab, ...)
  return(gt)
}

