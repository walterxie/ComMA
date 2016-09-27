# Extended bar chart 
# Author: Walter Xie
# Accessed on 15 May 2016

#' Percentage Bar
#' 
#' Percentage bar chart coloured by groups, which is extended from \code{\link{ggBarChart}}. 
#' 
#' @param df.to.melt A data frame required to \code{\link{melt}} before making a percent bar chart,
#' where Phyla have to be in row.names, otherwise it will be melted into values in y-axis to cause 
#' "Error: Discrete value supplied to continuous scale". 
#' For example, \url{https://github.com/walterxie/ComMA/blob/master/data-raw/reads.phyla.txt}.
#' @param melt.id A column name to \code{\link{melt}} and used to assign the colours.
#' @param title Graph title
#' @param x.lab,y.lab The label of x-axis or y-axis, such as plot names.
#' @param low, high Refer to \pkg{ggplot2} \code{\link{scale_fill_gradient}}. 
#' Default to low="white", high="steelblue".
#' @param autoWidth If TRUE, then use number of bars and legend columns 
#' to estimate pdf width automatically. Default to TRUE.
#' @param ... Other arguments passed to \code{\link{ggBarChart}}.
#' @keywords graph
#' @export
#' @examples 
#' data(reads.phyla)
#' bar.chart <- ggPercentageBarChart(reads.phyla, melt.id="TaxaGroup")
#' bar.chart$gg.plot
ggPercentageBarChart <- function(df.to.melt, melt.id, title="Percentage Bar Chart", 
                                 x.lab="", y.lab="", palette=NULL, x.meta.data=NULL,
                                 x.text.angle=90, autoWidth=TRUE, ...) {
  if (!is.element(tolower(melt.id), tolower(colnames(df.to.melt))))
    stop(paste0("Data frame column names do NOT have \"", melt.id, "\" for melt function !"))
  
  suppressMessages(require(reshape2))
  df.melt <- melt(df.to.melt, id=c(melt.id))
  #df.melt[,"variable"] <- factor(df.melt[,"variable"], levels = sort(unique(df.melt[,"variable"])))
  
  # sort legend
  legend.ord <- as.character(sort(unique(df.to.melt[,melt.id])))
  # move unclassified group to the last of legend 
  id.match <- grep("unclassified", legend.ord, ignore.case = TRUE)
  if (length(id.match) > 0)
    legend.ord <- legend.ord[c(setdiff(1:length(legend.ord), id.match),id.match)]
  df.melt[,melt.id] <- factor(df.melt[,melt.id], levels = legend.ord)
  
  if (! is.null(x.meta.data)) {
    id.match <- match(df.melt$variable, row.names(x.meta.data))
    if (length(id.match) < 1)
      stop("Invalid x.meta.data in row.names or 1st column !")
    df.melt$facet <- x.meta.data[,1]
    y.facet.id <- "facet"
  } else {
    y.facet.id <- NULL
  }
  
  if (! is.null(palette)) {
    pale <- palette
  } else {
    pale <- ComMA::getMyPalette(length(legend.ord))
    if (length(legend.ord) > length(pale)) {
      suppressMessages(require(colorspace))
      pale <- rainbow_hcl(length(legend.ord))
    }
  }
  # number of columns for legend
  legend.col = ceiling(length(legend.ord) / 25)
  
  p <- ComMA::ggBarChart(df.melt, x.id="variable", y.id="value", fill.id=melt.id, 
                         bar.pos="fill", y.trans="per", palette=pale, 
                         title=title, x.lab=x.lab, y.lab=y.lab, y.facet.id=y.facet.id,
                         x.text.angle=x.text.angle, legend.col=legend.col, ...)
  
  if (autoWidth)
    pdf.width = 1 + legend.col*2.5 + length(unique(df.melt[,"variable"])) * 0.2
  
  # Return a list containing the filename
  list(
    pdf.width = pdf.width,
    legend.len = length(legend.ord), # the length of legend
    legend.ord = legend.ord, # the legend in the order
    gg.plot = p # ggplot
  )
}

#' Group Abundance Bar
#' 
#' Group abundance bar chart for mulit-sample or mulit-dataset and coloured by groups, 
#' which is extended from \code{\link{ggBarChart}}. 
#' 
#' \strong{Note:} the coordinate is flipped in this chart. 
#' please make sure you recognise which is x or y,  
#' when defining the related arguments.   
#' 
#' @param df.to.melt A data frame required to \code{\link{melt}} before making the chart,
#' where Phyla have to be in row.names, otherwise it will be melted into values in y-axis to cause 
#' "Error: Discrete value supplied to continuous scale". But it  
#' For example, \url{https://github.com/walterxie/ComMA/blob/master/data-raw/reads.phyla.txt}.
#' @param melt.id A column name to \code{\link{melt}} and used to label axis.
#' @param colour.id The column name to give colours, and it will join 
#' \code{melt.id} as id to \code{\link{melt}}. Default to NULL.
#' @param prop.thre Make "Others" for any row total abundence < 
#' proportion threshold of the total of matrix abundence. 
#' Recommend 0.001 (0.1\%) if too many rows. Default to 0.
#' @param melt.levels,colour.levels The \code{\link{levels}} encode 
#' \code{\link{factor}} on melt.id and colour.id, which make data frame 
#' grouped for plotting. The default is a function to \code{\link{sort}}.
#' It can also use \code{\link{unique}} only for no sorting. 
#' Note: do not use \code{sort(,decreasing = T)}, it will ruin \code{factor}. 
#' Use \code{\link{rev}} instead.
#' @param title Graph title
#' @param x.lab,y.lab The label of x-axis or y-axis, such as plot names.
#' @param autoSize If TRUE, then use number of bars and legend columns 
#' to estimate pdf width and height automatically. Default to TRUE.
#' @param ... Other arguments passed to \code{\link{ggBarChart}}.
#' @keywords graph
#' @export
#' @examples 
#' data(reads.phyla)
#' reads.phyla$phylum <- rownames(reads.phyla)
#' bar.chart <- ggGroupAbundanceBar(reads.phyla, melt.id="phylum", colour.id="TaxaGroup", legend.row=1)
#' bar.chart$gg.plot
ggGroupAbundanceBar <- function(df.to.melt, melt.id, colour.id=NULL, prop.thre=0, 
                                melt.levels=function(x) rev(unique(x)), 
                                colour.levels=function(x) sort(unique(x)),
                                y.trans="log", legend.row=0, palette=NULL, 
                                title="Abundance Bar Chart", x.lab="", y.lab="Reads", 
                                autoSize=TRUE, verbose=TRUE, ...) {
  if (!is.element(tolower(melt.id), tolower(colnames(df.to.melt))))
    stop(paste0("Data frame column names do NOT have \"", melt.id, "\" for melt function !"))
  
  #####  extract abundence matrix #####
  melt.col.ind <- match(tolower(melt.id), tolower(colnames(df.to.melt)))
  abun.col.ind <- setdiff(1:ncol(df.to.melt), melt.col.ind)
  total.col.ind <- match("total", tolower(colnames(df.to.melt)))
  if (! is.na(total.col.ind))
    abun.col.ind <- setdiff(abun.col.ind, total.col.ind)
  
  if (! is.null(colour.id)) {
    colour.col.ind <- match(tolower(colour.id), tolower(colnames(df.to.melt)))
    if (is.na(colour.col.ind))
      stop("Invalid colour.id,", colour.id,  "not exsit in column names !\n")
    abun.col.ind <- setdiff(abun.col.ind, colour.col.ind)
  } 
  
  if (verbose)
    cat("The abundence data columns (exclude total) are ", 
        paste(colnames(df.to.melt)[abun.col.ind], collapse = ","), "\n")
  
  #####  make "Others" category #####
  # make "Others" for any row total abundence < proportion threshold of the total of matrix abundence
  if (prop.thre > 0) {
    if (prop.thre >= 1)
      stop("Invalid proportion threshold", prop.thre,  ", it must > 0 and < 1 !\n")
    total.thre <- sum(df.to.melt[ ,abun.col.ind]) * prop.thre
    if (verbose)
      cat("The total of matrix abundence = ", sum(df.to.melt[ ,abun.col.ind]), 
          ", proportion threshold = ", prop.thre, "\n")
    
    other.row.ind <- which(rowSums(df.to.melt[ ,abun.col.ind]) < total.thre)
    if (verbose)
      cat("Move", length(other.row.ind), "rows to 'Others'.\n")
    others <- colSums(df.to.melt[other.row.ind, abun.col.ind])
    if (length(others) != length(abun.col.ind))
      stop("'Others'", length(others),  "elements do not match", length(abun.col.ind), 
           "columns in abundence matrix !\n")
    
    df.to.melt <- df.to.melt[-other.row.ind,]
    # replicate last row
    df.to.melt <- rbind(df.to.melt, df.to.melt[nrow(df.to.melt),])
    # fill in "Others"
    rownames(df.to.melt)[nrow(df.to.melt)] <- "Others"
    df.to.melt[nrow(df.to.melt),melt.col.ind] <- "Others"
    df.to.melt[nrow(df.to.melt),abun.col.ind] <- others
    if (! is.na(total.col.ind)) {
      others.total <- sum(others)
      df.to.melt[nrow(df.to.melt),total.col.ind] <- others.total
    }
    if (! is.null(colour.id)) {
      df.to.melt[nrow(df.to.melt),colour.col.ind] <- "Others"
    }   
    
    if (nchar(x.lab) < 1)
      x.lab <- paste(nrow(df.to.melt)-1, "of", nrow(df.to.melt)-1+length(other.row.ind), 
                     "taxa at", melt.id[1], "or higher-level ( >", prop.thre*100, "% of total)")	
  }
  
  suppressMessages(require(reshape2))
  # colour.id == NULL still works
  df.melt <- melt(df.to.melt, id=c(melt.id, colour.id))
  #df.melt[,"value"] <- as.numeric(df.melt[,"value"])
  #df.melt[,"variable"] <- factor(df.melt[,"variable"], levels = sort(unique(df.melt[,"variable"])))
  
  # legend
  legend.ord <- c()
  if (! is.null(colour.id)) {
    # sort legend
    legend.ord <- colour.levels(df.melt[,colour.id])
    if (verbose)
      cat("Find", length(legend.ord), "unique groups for colouring.\n")
    
    df.melt[,colour.id] <- factor(df.melt[,colour.id], levels = legend.ord)
  }
  
  # default rev(df.melt[,melt.id]), note sort(,decreasing = T) will ruin factor 
  df.melt[,melt.id] <- factor(df.melt[,melt.id], levels = melt.levels(df.melt[,melt.id]))
  
  # number of columns for legend
  if (legend.row == 0)
    legend.row = ceiling(length(legend.ord) / 7)
  
  if (prop.thre == 0 && nchar(x.lab) < 1)
    x.lab <- paste(nrow(df.to.melt), "taxa at", melt.id[1], "or higher-level")	
  
  if (! is.null(palette)) {
    pale <- palette
  } else {
    pale <- ComMA::getMyPalette(length(legend.ord))
    if (length(legend.ord) > length(pale)) {
      suppressMessages(require(colorspace))
      pale <- rainbow_hcl(length(legend.ord))
    }
  }
  
  p <- ComMA::ggBarChart(df.melt, x.id=melt.id[1], y.id="value", fill.id=colour.id, 
                         y.facet.id="variable", bar.pos="identity", coord.flip=T,
                         y.trans=y.trans, title=title, x.lab=x.lab, y.lab=y.lab, 
                         palette=pale, legend.row=legend.row, verbose=verbose,
                         legend.position="top", legend.direction="horizontal", ...)
  
  n.row = nrow(df.to.melt)
  n.col = length(abun.col.ind)
  if (! is.na(total.col.ind))
    n.col = n.col + 1
  
  if (autoSize) {
    maxLabelLen= max(nchar( as.character(df.to.melt[,melt.id]) )) 
    pdf.width = 0.1 + maxLabelLen / 10 + (n.col-2) * 1.5
    pdf.height = 1 + legend.row * 1 + n.row * 0.12
  }
  
  # Return a list containing the filename
  list(
    pdf.width = pdf.width,
    pdf.height = pdf.height,
    legend.row = legend.row, # the length of legend
    legend.ord = legend.ord, # the legend in the order
    n.row = n.row, # number of taxa
    n.col = n.col, # include 'taxa' and 'group' col
    gg.plot = p # ggplot
  )
}


#' Y Across X Bar Chart
#' 
#' The bar chart shows the number of OTUs/reads in the bar, 
#' which simultaneously appeared at the number of samples in the value of x-axis. 
#' Given a data matrix whose structure is same as community matrix defined in \code{\link{ComMA}}, 
#' this function uses \code{\link{cmYAcrossX}} to aggregate it into another abundance matrix, 
#' and uses \code{\link{ggBarChart}} to plot.
#' 
#' @param community.matrix Community matrix (OTU table), where rows are 
#' OTUs or individual species and columns are sites or samples. 
#' @param terms The terms of matrix for x, y, and cells. 
#' Default to \code{c("OTUs", "samples", "reads")}.
#' @param title Graph title
#' @param x.lab,y.lab The label of x-axis or y-axis, such as plot names.
#' @param low, high Refer to \pkg{ggplot2} \code{\link{scale_fill_gradient}}. 
#' Default to low="white", high="steelblue".
#' @param autoWidth If TRUE, then use number of bars and legend columns 
#' to estimate pdf width automatically. Default to TRUE.
#' @param ... Other arguments passed to \code{\link{ggBarChart}}.
#' @keywords graph
#' @export
#' @examples  
#' community.matrix <- getCommunityMatrix("16S", isPlot=TRUE, minAbund=1)
#' bar.yx <- ggYAcrossXBar(community.matrix)
ggYAcrossXBar <- function(community.matrix, terms=c("OTUs", "samples", "reads"),
                          title="The number of OTUs/reads across the number of samples", 
                          x.lab="Number of samples crossed", y.lab="Number of OTUs/reads",
                          y.trans="log", auto.scale.y=TRUE, x.scale="discrete", x.interval=1, 
                          x.text.angle=0, legend.title="", ...) {
  cm.aggre <- ComMA::cmYAcrossX(community.matrix, terms=terms)
  suppressMessages(suppressWarnings(require(reshape2)))
  y.term <- terms[2]
  df <- melt(cm.aggre, id=y.term)
  
  if (x.scale=="discrete") {
    df[,y.term] <- as.character(df[,y.term])
    df[,y.term] <- factor(df[,y.term], unique(df[,y.term]))
  }
  
  p <- ComMA::ggBarChart(df, x.id=y.term, y.id="value", fill.id="variable", 
                         y.trans=y.trans, auto.scale.y=auto.scale.y, 
                         title=title, x.lab=x.lab, y.lab=y.lab, 
                         x.text.angle=x.text.angle, x.scale=x.scale, x.interval=x.interval, 
                         legend.title=legend.title, ...)
  
  return(p)
}

