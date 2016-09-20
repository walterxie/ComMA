# Author: Andrew Dopheide, Walter Xie
# Accessed on 16 Nov 2016


### Get OTU and sequence counts by taxonomic group ###

#' @name countOTUsReads
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
#' \code{getCountSums} return a list of counts of OTUs and reads given a list of
#' data sets.
#' 
#' @param ... Input of a list of community matrices, or comma separated multi-inputs.
#' @param input.list Default to TRUE to unwrap list(...) to 
#' get the actual list if the input is a list of cm. 
#' @param taxa.rank The rank to display in the axis for counts of OTUs and reads.
#' @param group.rank The rank to group \code{taxa.rank} to colour the figure. 
#' @keywords taxonomy
#' @export
#' @examples 
#' all.counts.sums <- ComMA::getCountSums(cm.taxa.list, input.list=T, group.rank="kingdom", taxa.rank="phylum")
#' 
#' @rdname countOTUsReads
getCountSums <- function(..., input.list=FALSE, group.rank="kingdom", taxa.rank="phylum"){
  cm.taxa.list <- validateInputList(..., input.list=input.list) 
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
    require(data.table)
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
  require(data.table)
  return(rbindlist(counts.sums.list)) # list
}


#' @details 
#' \code{plotTaxonomy} makes a plot of OTU and sequence counts by taxonomic group 
#' sumarised by \code{getCountSums}.
#' 
#' Diamonds represent the number of sequences, 
#' open circles the number of OTUs including singleton OTUs, 
#' and filled circles the number of OTUs excluding singleton OTUs.
#' 
#' @param gene.level The level to order 'gene' column.  
#' @param group.level The level to order 'gene' column.
#' @param x.lab,y.lab Label for x or y axis.
#' @param taxa.ref The taxonomic reference data set to order taxonomic names nicely. 
#' Default to "" to ignore it, the logic is length(taxa.ref) > 0.
#' @keywords taxonomy
#' @export
#' @examples 
#' plotTaxonomy(all.counts.sums, taxa.ref=taxa.ref)
#' 
#' @rdname countOTUsReads
plotTaxonomy <- function(all.counts.sums, taxa.ref="", taxa.rank="phylum", group.rank="kingdom", 
                         gene.level=c("16S", "18S", "26S", "ITS", "COI-300", "COI-650"),
                         group.level=c("ARCHAEA","BACTERIA","EUKARYOTA","PROTOZOA","CHROMISTA","FUNGI","PLANTAE","ANIMALIA","Unknown"),
                         x.lab="Phylum (or higher-level taxon)", y.lab="Number of sequences or OTUs"){
  # z is a list
  z <- all.counts.sums[all.counts.sums[,taxa.rank] != 0, ]
  # fix name for convience
  colnames(z)[match(c(taxa.rank, group.rank), colnames(z))] <- c("taxa", "group")
  if (length(taxa.ref) > 0) {
    # make case insensitive
    colnames(taxa.ref) <- tolower(colnames(taxa.ref))
    if (! taxa.rank %in% colnames(taxa.ref))
      stop("Invalid taxonomic reference data set: no ", taxa.rank, " column !")
    z$taxa <- taxa.ref[match(tolower(z$taxa), tolower(taxa.ref[,taxa.rank])), taxa.rank]
  }
  require(reshape2)
  z <- as.data.frame(z)
  z <- melt(z, id.vars = c("taxa", "group", "gene"))
  # Order factors
  if (length(taxa.ref) > 0)
    z$taxa <- factor(z$taxa, ordered = TRUE, levels = rev(unique(taxa.ref[,taxa.rank])))
  z$gene <- factor(z$gene, ordered = TRUE, levels = gene.level)
  z$group <- gsub("root|No hits|Not assigned|cellular organisms", "Unknown", z$group)
  z$group <- factor(z$group, ordered = TRUE, levels = group.level)
  z <- na.omit(z)
  #print(z)
  require(ggplot2)
  p <- ggplot(z) + 
    geom_point(data = z[z$variable == "count.1",], aes(x = taxa, y = value, colour = group), shape = 1, size = 2) +
    geom_point(data = z[z$variable == "sum.1",], aes(x = taxa, y = value, colour = group), shape = 5, size = 2, show.legend = NA) +
    geom_point(data = z[z$variable == "count.2",], aes(x = taxa, y = value, colour = group), shape = 16, size = 2) +
    #geom_point(data = z[z$variable == "sums.2",], aes(x = taxa, y = value, colour = group), shape = 17, size = 2) +
    geom_line(data = z, aes(x = taxa, y = value, group=interaction(taxa, gene), colour = group), size = 0.5, alpha = 0.5) + 
    facet_grid( ~ gene) + guides(fill = guide_legend(reverse = FALSE)) +
    coord_flip() + theme_bw() + ylab(y.lab) + xlab(x.lab) + 
    scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000), label=ComMA::scientific_10) +
    #scale_color_manual(values = pal(length(unique(z$group)))) +
    theme(strip.background = element_blank(), plot.title = element_text(size = 9),
          #plot.margin=unit(c(0.2,0.5,0.2,0.8), "cm"), panel.margin = unit(0.8, "lines"), 
          axis.title.y=element_text(vjust=2), axis.title.x=element_text(vjust=-0.2), 
          axis.text.x = element_text(vjust=-0.1)) 
  return(p)
}

#' @details 
#' \code{sumReadsOTUs} is a high-level function to use \code{getCountSums} sumarise counts,
#' and call \code{plotTaxonomy} to plot figure.
#' 
#' @param cm.taxa.list A list of \code{cm.taxa} merged by \code{\link{mergeCMTaxa}}. 
#' Require proper names of the list.  
#' @param col.ranks A vector or string of column name(s) of taxonomic ranks in the taxa table, 
#' detail to see \code{\link{mergeCMTaxa}}. 
#' Default to c("kingdom", "phylum", "class", "order").
#' @param fig.folder Folder to contain the figure in PDF from \code{plotTaxonomy}.
#' @param pdf.width,pdf.height,units Parameters in \code{\link{ggsave}} to adjust figure.
#' @param table.folder Folder to save the result from \code{getCountSums}.
#' @keywords taxonomy
#' @export
#' @examples 
#' 
#' @rdname countOTUsReads
sumReadsOTUs <- function(cm.taxa.list, taxa.ref="", taxa.rank="phylum", group.rank="kingdom", 
                         col.ranks=c("kingdom", "phylum", "class", "order"), 
                         gene.level=c("16S", "18S", "26S", "ITS", "COI-300", "COI-650"),
                         group.level=c("ARCHAEA","BACTERIA","EUKARYOTA","PROTOZOA","CHROMISTA","FUNGI","PLANTAE","ANIMALIA","Unknown"),
                         x.lab="Phylum (or higher-level taxon)", y.lab="Number of sequences or OTUs",
                         fig.folder="./figures", table.folder="./outputs",
                         pdf.width = 260, pdf.height = 200, units = "mm"){
  all.counts.sums <- ComMA::getCountSums(cm.taxa.list, input.list=T, taxa.rank=taxa.rank, group.rank=group.rank)
  
  if (!is.na(table.folder)) {
    write.table(all.counts.sums, file = file.path(table.folder, paste0("Overall_counts_sums_by_", taxa.rank, ".txt")), 
                sep = "\t", quote = FALSE, col.names = NA)
  }
  if (!is.na(fig.folder)) {
    p <- ComMA::plotTaxonomy(all.counts.sums, taxa.ref=taxa.ref, gene.level=gene.level, 
                             group.level=group.level, x.lab=x.lab, y.lab=y.lab)
    ggsave(p, file = file.path(fig.folder, paste0("Overall_taxonomy_OTUs_reads_by_", taxa.rank, ".pdf")), 
           width = pdf.width, height = pdf.height, units = units)
  }
  
  return(all.counts.sums) 
}
