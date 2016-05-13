# Extended bar chart 
# Author: Walter Xie
# Accessed on 15 May 2016

#' Percentage bar chart coloured by groups, which is extended from \code{\link{ggBarChart}}. 
#' 
#' @param df.to.melt A data frame required to \code{\link{melt}} before making a percent bar chart. 
#' For example,
#' \tabular{rrrrr}{
#'   Phyla \tab 16s \tab 18s \tab ITS \tab TaxaGroup\cr
#'   Actinobacteria \tab 958 \tab 1 \tab 3 \tab Bacteria\cr
#'   Crenarchaeota \tab 1 \tab 0 \tab 0 \tab Archaea\cr
#'   Ascomycota \tab 2 \tab 765 \tab 971 \tab Fungi 
#' } 
#' @param melt.id A column name to \code{\link{melt}} and used to assign the colours.
#' @param title Graph title
#' @param x.lab,y.lab The label of x-axis or y-axis, such as plot names.
#' @param low, high Refer to \pkg{ggplot2} \code{\link{scale_fill_gradient}}. 
#' Default to low="white", high="steelblue".
#' @param autoWidth If TRUE, then use number of bars and legend columns 
#' to estimate pdf width automatically. Default to TRUE.
#' @param ... More from \code{\link{ggBarChart}}.
#' @keywords graph
#' @export
#' @examples 
#' data(reads.phyla)
#' reads.phyla
#' bar.chart <- ggPercentageBarChart(reads.phyla, melt.id="TaxaGroup")
#' bar.chart$gg.plot
#' pdf.ggplot(bar.chart$gg.plot, fig.path="taxa-percentage-bar.pdf", width=bar.chart$pdf.width, height=8) 
ggPercentageBarChart <- function(df.to.melt, melt.id, title="Percentage Bar Chart", 
                                 x.lab="", y.lab="", palette=NULL, 
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
                         bar.pos="fill", palette=pale, 
                         x.text.angle=x.text.angle, y.trans="per", 
                         title=title, x.lab=x.lab, y.lab=y.lab, 
                         legend.col=legend.col, ...)
  
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


ggAbundanceBar <- function(df.to.melt, melt.id, title="Abundance Bar Chart", 
                                 x.lab="", y.lab="", palette=NULL, 
                                 x.text.angle=90, autoWidth=TRUE, ...) {
  if (!is.element(tolower(melt.id), tolower(colnames(df.to.melt))))
    stop(paste0("Data frame column names do NOT have \"", melt.id, "\" for melt function !"))
  
  suppressMessages(require(reshape2))
  df.melt <- melt(df.to.melt, id=c(melt.id))
  #df.melt[,"variable"] <- factor(df.melt[,"variable"], levels = sort(unique(df.melt[,"variable"])))
  
  #####  make "Others" category #####
  # make "Others" by percThr: percentage threshold of total reads in taxaAssg[ ,total_string]
  if (percThr > 0) {
    totalThr <- sum(taxaAssg[ ,total_string]) * percThr
    Others <- colSums(taxaAssg[which(taxaAssg[ ,total_string]<=totalThr),2:colTotal])
    Others <- c("Others",Others,rep("Others", ncol(taxaAssg)-colTotal),"")	
    taxaAssgId <- which(taxaAssg[ ,total_string]>totalThr)
    
    cat("Filter out", nrow(taxaAssg)-length(taxaAssgId), "taxa to Others category whose total <=", 
        percThr, "of the total of whole matrix.\n")  
    
    xlab <- xlab( paste(length(taxaAssgId), " of ", nrow(taxaAssg), rankLevel, " (", tolower(y_string), 
                        " > ", percThr*100, "% of total) ", sep = "") )		
    
    # avoid error: invalid factor level, NA generated
    taxaAssg <- data.frame(lapply(taxaAssg, as.character), check.names=FALSE, stringsAsFactors=FALSE)
    taxaAssg <- rbind(taxaAssg[taxaAssgId,], Others)
    
    group = unique(taxaAssg[,groupLevel]) # include Other
    # -1 to consider "Others" category in the last row
    cat("Select", nrow(taxaAssg)-1, rankLevel, "by percThr = ", percThr, ", assigned to", length(group), groupLevel, 
        ", total", y_string, "=",  sum(taxaAssg[ ,total_string]), ".\n")  
  }
  
  
  
  
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


#' The bar chart shows the number of OTUs in the bar across the number of samples in the value of x-axis, 
#' which is extended from \code{\link{ggBarChart}}. 
#' 
#' @param community.matrix Community matrix (OTU table), where rows are 
#' OTUs or individual species and columns are sites or samples. See \code{\link{ComMA}}.
#' @param title Graph title
#' @param x.lab,y.lab The label of x-axis or y-axis, such as plot names.
#' @param low, high Refer to \pkg{ggplot2} \code{\link{scale_fill_gradient}}. 
#' Default to low="white", high="steelblue".
#' @param autoWidth If TRUE, then use number of bars and legend columns 
#' to estimate pdf width automatically. Default to TRUE.
#' @param ... More from \code{\link{ggBarChart}}.
#' @keywords graph
#' @export
#' @examples  
#' community.matrix <- getCommunityMatrix("16S", isPlot=TRUE, minAbund=1)
#' bar.yx <- ComMA::ggBarYAcrossX(community.matrix)
ggBarYAcrossX <- function(community.matrix, title="The number of OTUs/reads across the number of samples", 
                          x.lab="Number of samples crossed", y.lab="Number of OTUs/reads",
                          y.trans="log", auto.scale.y=TRUE, x.scale="discrete", x.interval=1, 
                          x.text.angle=0, legend.title="", ...) {
  cm.aggre <- ComMA::cmYAcrossX(community.matrix)
  suppressMessages(require(reshape2))
  df <- melt(cm.aggre, id="samples")
  
  if (x.scale=="discrete") {
    df[,"samples"] <- as.character(df[,"samples"])
    df[,"samples"] <- factor(df[,"samples"], unique(df[,"samples"]))
  }
  
  p <- ComMA::ggBarChart(df, x.id="samples", y.id="value", fill.id="variable", 
                         y.trans=y.trans, auto.scale.y=auto.scale.y, 
                         title=title, x.lab=x.lab, y.lab=y.lab, 
                         x.text.angle=x.text.angle, x.scale=x.scale, x.interval=x.interval, 
                         legend.title=legend.title, ...)
  
  return(p)
}

