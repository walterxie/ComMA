# Author: Andrew Dopheide, Walter Xie
# Accessed on 16 Nov 2016


### Get OTU and sequence counts by taxonomic group ###

#' @name ReadsOTUsTaxa
#' @title Taxonomy analyses for multiple data sets
#'
#' @description Taxonomic summary of OTUs and reads for multiple data sets, 
#' and bar chart for visualization.
#' "count.1" column is number of OTUs including singletons, 
#' "sum.1" is number of reads including singletons, 
#' "count.2" is OTUs excluding singletons, 
#' "sum.2" is reads excluding singletons. 
#' 
#' @details 
#' \code{getCountsSums} return a data frame
#' concatenating counts of OTUs and reads given a list of data sets.
#' 
#' @param ... Input of a list of community matrices, or comma separated multi-inputs.
#' @param input.list Default to TRUE to unwrap list(...) to 
#' get the actual list if the input is a list of cm. 
#' @param taxa.rank The rank to display in the axis for counts of OTUs and reads.
#' @param group.rank The rank to group \code{taxa.rank} to colour the figure. 
#' @keywords taxonomy
#' @export
#' @examples 
#' all.counts.sums <- ComMA::getCountsSums(cm.taxa.list, input.list=T, group.rank="kingdom", taxa.rank="phylum")
#' 
#' @rdname ReadsOTUsTaxa
getCountsSums <- function(..., input.list=FALSE, taxa.rank="phylum", group.rank="kingdom"){
  cm.taxa.list <- unwrapInputList(..., input.list=input.list) 
  cat("Count OTU/Reads on", length(cm.taxa.list), "data sets.\n") 
  
  counts.sums.list <- list()
  for(n in 1:length(cm.taxa.list)){
    cm.taxa <- cm.taxa.list[[n]]
    attr.cm.ta <- attributes(cm.taxa)
    ncol.cm <- attr.cm.ta$ncol.cm
    col.ranks <- attr.cm.ta$col.ranks
    
    if (is.null(ncol.cm) || is.null(col.ranks)) 
      stop("Input cm.taxa", names(cm.taxa.list)[n], "should have attributes ncol.cm and col.ranks !")
    
    if (ncol.cm > 1) {
      total <- rowSums(cm.taxa[,1:ncol.cm]) # do not have 'total' column
    } else {
      total <- cm.taxa[,1] # only have 'total' column
    }
    OTU.sums <- cbind(cm.taxa[,(ncol.cm+1):ncol(cm.taxa)], total)
    colnames(OTU.sums)[ncol(OTU.sums)] <- "total"
    
    sums.1 <- aggregate(as.formula(paste("total ~", taxa.rank, "+", group.rank)), 
                        data=OTU.sums[OTU.sums$total > 0,], FUN=sum)
    sums.2 <- aggregate(as.formula(paste("total ~", taxa.rank, "+", group.rank)), 
                        data=OTU.sums[OTU.sums$total > 1,], FUN=sum)
    counts.1 <- aggregate(as.formula(paste("total ~", taxa.rank, "+", group.rank)), 
                          data=OTU.sums[OTU.sums$total > 0,], FUN=function(x) sum(x>0))
    counts.2 <- aggregate(as.formula(paste("total ~", taxa.rank, "+", group.rank)), 
                          data=OTU.sums[OTU.sums$total > 1,], FUN=function(x) sum(x>0))
    
    counts.sums <- merge(counts.1, sums.1, by = c(taxa.rank, group.rank))
    counts.sums.2 <- merge(counts.2, sums.2, by = c(taxa.rank, group.rank))
    counts.sums <- merge(counts.sums, counts.sums.2, by = c(taxa.rank, group.rank), all = TRUE)
    counts.sums[is.na(counts.sums)] <- 0
    colnames(counts.sums) <- c(taxa.rank, group.rank, "count.1", "sum.1", "count.2", "sum.2")
    counts.sums$gene <- names(cm.taxa.list)[n]
    
    label <- paste(names(cm.taxa.list)[n], taxa.rank, sep = ".")
    cat("Data set", label, "has", nrow(counts.sums), taxa.rank, 
        "and", length(unique(counts.sums[,group.rank])), group.rank,"\n")
    counts.sums.list[[label]] <- counts.sums
  }
  counts.sums.df <- data.frame(do.call("rbind", counts.sums.list))
  return(counts.sums.df)
}


#' @details 
#' \code{plotTaxonomy} plots a circle-diamond graph of reads and OTUs coloured by taxonomic group 
#' and sumarised by \code{getCountsSums}.
#' 
#' Diamonds represent the number of sequences, 
#' open circles the number of OTUs including singleton OTUs, 
#' and filled circles the number of OTUs excluding singleton OTUs.
#' 
#' @param all.counts.sums The output produced by \code{getCountsSums}.
#' @param taxa.ref The taxonomic reference data set to order taxonomic names nicely.
#' But it removes the taxonomy not exsiting in the reference.
#' Default to use a taxonomy reference data set \code{taxa.ref.PLOSONE.2015} in the package.
#' Set \code{taxa.ref=""}, where length(taxa.ref) > 1, to plot all taxonomy.
#' @param count.unclassified Count the taxonomy having 'unclassified' prefix, 
#' mostly in RDP classification, when using \code{taxa.ref.PLOSONE.2015}. The default is FALSE.
#' @param gene.levels The level to order 'gene' column.  
#' @param group.levels The level to order 'gene' column.
#' @param exclude.singletons If 0, then do not plot the number of 
#' sequences and OTUs excluding singleton OTUs. 
#' If 1, as default, just do not plot sequences excluding singleton. 
#' If 2, then plot all (including/excluding singleton).
#' @param x.lab,y.lab,title,title.size,legend.title,palette See \code{\link{ggPlot}}.
#' @keywords taxonomy
#' @export
#' @examples 
#' plotTaxonomy(all.counts.sums, taxa.ref=taxa.ref)
#' 
#' @rdname ReadsOTUsTaxa
plotTaxonomy <- function(all.counts.sums, taxa.ref=ComMA::taxa.ref.PLOSONE.2015, exclude.singletons=1,
                         taxa.rank="phylum", group.rank="kingdom", count.unclassified=FALSE,
                         gene.levels=c("16S", "18S", "26S", "ITS", "COI-300", "COI-650"),
                         group.levels=c("ARCHAEA","BACTERIA","EUKARYOTA","PROTOZOA","CHROMISTA","FUNGI","PLANTAE","ANIMALIA","Unknown"),
                         x.lab="Phylum (or higher-level taxon)", y.lab="Number of sequences or OTUs", 
                         legend.title=NULL, palette=NULL, title="", title.size = 10, verbose=TRUE){
  if (!taxa.rank %in% colnames(all.counts.sums) || !group.rank %in% colnames(all.counts.sums))
    stop("Cannot find 'taxa.rank' ", taxa.rank, " or 'group.rank'", group.rank,
         " in the input 'all.counts.sums' columns !")
  
  # in case if z is a data.table object from rbindlist
  z <- as.data.frame(all.counts.sums)
  # fix name for convience
  colnames(z)[match(c(taxa.rank, group.rank), colnames(z))] <- c("taxa", "group")
  if (length(taxa.ref) > 1) {
    # make case insensitive
    colnames(taxa.ref) <- tolower(colnames(taxa.ref))
    if (! taxa.rank %in% colnames(taxa.ref))
      stop("Invalid taxonomic reference data set: no ", taxa.rank, " column !")
    if (count.unclassified) {
      # change rdp unclassified to Not assigned to fit in factor
      z$taxa <- gsub("^unclassified$", "Not assigned", z$taxa, ignore.case = T)
      z$taxa <- gsub("^unclassified ", "", z$taxa, ignore.case = T)
    }
      
    z$taxa <- taxa.ref[match(tolower(z$taxa), tolower(taxa.ref[,taxa.rank])), taxa.rank]
  } 
  
  require(reshape2)
  z <- melt(z, id.vars = c("taxa", "group", "gene"))
  z$group <- gsub("root|No hits|Not assigned|cellular organisms|unclassified", "Unknown", z$group, ignore.case=T)
  # Order factors
  z$gene <- factor(z$gene, ordered = TRUE, levels = gene.levels)
  z$group <- factor(z$group, ordered = TRUE, levels = group.levels)
  if (length(taxa.ref) > 1) {
    taxa.levels <- rev(unique(taxa.ref[,taxa.rank]))
    if (count.unclassified) {
      # get unclassified key word back
      z$taxa <- gsub("Not assigned", "unclassified", z$taxa, ignore.case = T)
      taxa.levels[taxa.levels=="Not assigned"] <- "unclassified"
    }
    z$taxa <- factor(z$taxa, ordered = TRUE, levels = taxa.levels)
  } else {
    taxa.levels <- rev(unique(z[order(z$gene,z$group,z$taxa), "taxa"]))#rev(unique(z$taxa))
    z$taxa <- factor(z$taxa, levels=taxa.levels, ordered=TRUE)
  }
  
  z <- na.omit(z)
  #print(z)
  
  if (is.null(legend.title))
    legend.title <- group.rank
  require(ggplot2)
  p <- ggplot(z) +  
    geom_point(data = z[z$variable == "count.1",], aes(x = taxa, y = value, colour = group), shape = 1, size = 2) +
    geom_point(data = z[z$variable == "sum.1",], aes(x = taxa, y = value, colour = group), shape = 5, size = 2, show.legend = NA)
  if (exclude.singletons > 0)
    p <- p + geom_point(data = z[z$variable == "count.2",], aes(x = taxa, y = value, colour = group), shape = 16, size = 2)
  if (exclude.singletons > 1)
    p <- p + geom_point(data = z[z$variable == "sums.2",], aes(x = taxa, y = value, colour = group), shape = 17, size = 2)
  
  p <- p + geom_line(data = z, aes(x = taxa, y = value, group=interaction(taxa, gene), colour = group), size = 0.5, alpha = 0.5) +
    facet_grid( ~ gene) + guides(fill = guide_legend(reverse = FALSE)) +
    coord_flip() + theme_bw() + ylab(y.lab) + xlab(x.lab) + labs(colour=legend.title) + ggtitle(title) +
    scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000), label=ComMA::scientific_10)
  
  #scale_color_manual(values = pal(length(unique(z$group)))) +
  p <- ggOptPalette(p, palette=palette, verbose=verbose)
  
  p <- p + theme(strip.background = element_blank(), plot.title = element_text(size = title.size),
                 #plot.margin=unit(c(0.2,0.5,0.2,0.8), "cm"), panel.margin = unit(0.8, "lines"), 
                 axis.title.y=element_text(vjust=2), axis.title.x=element_text(vjust=-0.2), 
                 axis.text.x = element_text(vjust=-0.1)) 
  return(p)
}

#' @details 
#' \code{summReadsOTUsPipeline} is a pipeline to summarise reads and OTUs by \code{getCountsSums},
#' and plot a circle-diamond graph by \code{plotTaxonomy}.
#' 
#' @param cm.taxa.list A list of \code{cm.taxa} merged by \code{\link{mergeCMTaxa}}. 
#' Require proper names of the list.  
#' @param col.ranks A vector or string of column name(s) of taxonomic ranks in the taxa table, 
#' detail to see \code{\link{mergeCMTaxa}}. 
#' Default to c("kingdom", "phylum", "class", "order").
#' @keywords taxonomy
#' @export
#' @examples 
#' 
#' @rdname ReadsOTUsTaxa
summReadsOTUsPipeline <- function(cm.taxa.list, taxa.ref=ComMA::taxa.ref.PLOSONE.2015, taxa.rank="phylum", group.rank="kingdom", 
                                  col.ranks=c("kingdom", "phylum", "class", "order"), 
                                  gene.levels=c("16S", "18S", "26S", "ITS", "COI-300", "COI-650"),
                                  group.levels=c("ARCHAEA","BACTERIA","EUKARYOTA","PROTOZOA","CHROMISTA","FUNGI","PLANTAE","ANIMALIA","Unknown"),
                                  ...){
  all.counts.sums <- ComMA::getCountsSums(cm.taxa.list, input.list=T, taxa.rank=taxa.rank, group.rank=group.rank)
  
  p <- ComMA::plotTaxonomy(all.counts.sums, taxa.ref=taxa.ref, gene.levels=gene.levels, 
                           group.levels=group.levels, ...)
  
  list(ggplot=p, all.counts.sums=all.counts.sums) 
}
