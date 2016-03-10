# Graph
# Author: Walter Xie
# Accessed on 10 Mar 2016

#' Create a heat map using ggplot. 
#' 
#' @param df A data frame to \code{\link{melt}} and then make a heat map. 
#' For example,
#' \tabular{rrrr}{
#'   plot \tab 16s \tab 18s \tab ITS\cr
#'   CM30c39 \tab 2 \tab 1 \tab 3\cr
#'   CM30c44 \tab 10 \tab 26 \tab 15\cr
#'   Plot01 \tab 6 \tab 5 \tab 6 
#' } 
#' @param id.melt A column name as a \code{\link{factor}}.
#' @param fname.path The full path of image file.
#' @param title Graph title
#' @param x.lab, y.lab The label of x-axis or y-axis, such as plot names
#' @param low, high Refer to \pkg{ggplot2} \code{\link{scale_fill_gradient}}. Default to low="white", high="steelblue".
#' @param width, height Refer to \code{\link{pdf}}. Default to width=6, height=6.
#' @keywords graph
#' @export
#' @examples 
#' ranks.by.gourp <- data.frame(plot=c("Plot03","Plot02","Plot01"), `16s`=c(3,2,1), `18s`=c(1,2,3), ITS=c(2,1,3), check.names = F)
#' ranks.by.gourp
#' heatmapGgplot(df=ranks.by.gourp, id.melt="plot", fname.path="plot-prior-example-heatmap.pdf")
heatmapGgplot <- function(df, id.melt, fname.path, title="Heatmap", x.lab="", y.lab="", 
                             low="white", high="steelblue", width=6, height=6) {
  if (!is.element(id.melt, tolower(colnames(df))))
    stop(paste0("Data frame column names do NOT have \"", id.melt, "\" for melt function !"))
  
  require(reshape2)
  
  breaks.rank <- round(seq(1, nrow(df), length.out = 5), digits = 0)
  ranks.melt <- melt(df, id=c(id.melt))
  ranks.melt[,id.melt] <- factor(ranks.melt[,id.melt], levels=unique(ranks.melt[,id.melt]))
  
  require(ggplot2)
  # variable is all group names, such as "16S" or "FUNGI"
  # value is ranks for each group
  p <- ggplot(ranks.melt, aes_string(x="variable", y=id.melt)) + geom_tile(aes(fill=value)) + 
    scale_fill_gradient(na.value="transparent", low=low, high=high, name="rank", breaks=breaks.rank) +
    xlab(x.lab) + ylab(y.lab) + ggtitle(title) +
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
  
  pdf(fname.path, width=width, height=height)	
  print(p)
  invisible(dev.off()) 
}




