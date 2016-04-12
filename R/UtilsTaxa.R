# Author: Walter Xie
# Accessed on 11 Apr 2016

#' @name utilsTaxa
#' @title Utils to preprocess taxa table
#'
#' @description Utils to preprocess taxa table, 
#' and make it easy for visualization.
#' 
#' @details 
#' \code{subsetTaxaTable} takes or excludes a subset of given a taxa table at given rank.
#' 
#' @param taxa.table A data frame to contain taxonomic classifications of OTUs. 
#' Columns are taxonomy at the rank or lineage, rows are OTUs which need to 
#' match rows from community matrix. Use \code{\link{readTaxaTable}} to get it from file.
#' @param taxa.group The taxonomic group, the values can be 'all', 'assigned', or 
#' Group 'all' includes everything.
#' Group 'assigned' removes all uncertain classifications including 
#' 'root', 'cellular organisms', 'No hits', 'Not assigned'. 
#' Alternatively, any high-ranking taxonomy in your taxonomy file 
#' can be used as a group or multi-groups (seperated by "|"), 
#' such as 'BACTERIA', 'Proteobacteria', etc. But they have to be 
#' in the same rank column in the file. Default to remove all 
#' uncertain classifications, even when group(s) assigned.
#' @param rank The rank column in the file to specify where to 
#' search for \code{taxa.group}. 
#' @param include Define whether include or exclude given \code{taxa.group}. 
#' Default to TRUE.
#' @param verbose More details. Default to TRUE.
#' @keywords utils
#' @export
#' @examples 
#' tt.sub <- subsetTaxaTable(tt.megan, taxa.group="Proteobacteria", rank="phylum")
#' tt.sub <- subsetTaxaTable(tt.megan, taxa.group="Cnidaria|Brachiopoda|Echinodermata|Porifera", rank="phylum", include=FALSE)
#' 
#' @rdname utilsTaxa
subsetTaxaTable <- function(taxa.table, taxa.group="assigned", rank="kingdom", include=TRUE, verbose=TRUE) {
  if (include) {
    # include PROTISTS, taxa.group="CHROMISTA|PROTOZOA", rank="kingdom"
    taxa.table <- subset(taxa.table, (grepl(taxa.group, taxa.table[,rank], ignore.case = T))) 
  } else { 
    # exclude some phyla, taxa.group="Cnidaria|Brachiopoda|Echinodermata|Porifera", rank="phylum"
    taxa.table <- subset(taxa.table, !grepl(taxa.group, taxa.table[,rank], ignore.case = T)) 
  }
  return(taxa.table)
}


#' @details 
#' \code{assignTaxa} provides a list of taxonomic assignments with abudence 
#' from community matrix. The function is iterated through \code{col.names}, 
#' and \code{\link{aggregate}}s abudence into taxonomy based on the rank in 
#' \code{col.names}. All sequences either classified as 
#' "root|cellular organisms|No hits|Not assigned|unclassified sequences"
#' from BLAST + MEGAN, or confidence < \emph{min.conf} threshold from RDP, 
#' are changed to "unclassified", which will be moved to the last row.  
#' 
#' @param community.matrix Community matrix (OTU table), where rows are 
#' OTUs or individual species and columns are sites or samples. See \code{\link{ComMA}}.
#' @param classifier The classifier is used to generate \code{taxa.table}. 
#' Value is MEGAN or RDP. Default to MEGAN.
#' @param min.conf The confidence threshold to drop rows < \emph{min.conf}.
#' @param has.total If 0, then only return abudence by samples (columns) of community matrix. 
#' If 1, then only return toal abudence. If 2, then return abudence by samples (columns) and total. 
#' Default to 1.
#' @param col.names A vector or string of column name(s) of taxonomic ranks in the taxa table, 
#' which will determine the aggregated abundence matrix. They have to be full set or subset of 
#' c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"). 
#' Default to c("kingdom", "phylum", "class", "order", "family", "genus").
#' @keywords utils
#' @export
#' @examples 
#' ta.megan <- assignTaxa(community.matrix, tt.megan)
#' ta.rdp <- assignTaxa(community.matrix, tt.rdp, col.names="kingdom", classifier="RDP")
#' 
#' @rdname utilsTaxa
assignTaxa <- function(community.matrix, taxa.table, classifier="MEGAN", min.conf=0.8, has.total=1,
                       col.names=c("kingdom", "phylum", "class", "order", "family", "genus")) {
  ranks <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  if (length(col.names) < 1 || !all(col.names %in% ranks)) 
    stop(cat("Invaild column names for ranks !\nUse one element or a subset of", 
             paste(ranks, collapse = ",")))
  if (!all(col.names %in% colnames(taxa.table))) 
    stop(cat("Column names in taxa.table do not have", 
             paste(col.names, collapse = ",")))
  if ("confidence" %in% colnames(taxa.table) && classifier != "RDP")
    stop("Find 'confidence' column, please define classifier = 'RDP' !")
  
  ###### merge community.matrix, taxa.table 
  cm <- data.frame(row.names = rownames(community.matrix))
  if (has.total!=1) 
    cm <- community.matrix
  if (has.total > 0) 
    cm$Total <- rowSums(community.matrix) 
  ncol.cm <- ncol(cm)
  
  if (classifier == "RDP") {
    if (! "confidence" %in% colnames(taxa.table))
      stop("Invaild file format from RDP: column 'confidence' is required !")
    
    ### drop rows < min.conf
    n.row.un <- nrow(taxa.table[taxa.table[,"confidence"] < min.conf, ])
    taxa.table[taxa.table[,"confidence"] < min.conf, col.names] <- "unclassified"
    
    cat("Set", n.row.un, "rows as unclassified from the total of", nrow(taxa.table), 
        "in RDP taxa table, whose confidence <", min.conf, ".\n")
  }
  
  # taxa.assign 1st col is "row.names", "ncol.cm" columns abundence, and length(col.names) columns rank
  taxa.assign <- merge(cm, taxa.table, by = "row.names")
  # order rank by rank
  ord.cmd = parse(text = paste0('taxa.assign[order(taxa.assign[,"', paste(col.names, collapse = '"], taxa.assign[,"'), '"]),]')) 
  taxa.assign <- eval(ord.cmd)
  
  cat("Merge", nrow(cm), "rows in community matrix with", nrow(taxa.table), "rows in taxa table, get", 
      nrow(taxa.assign), "classifications.\n")
  
  ###### aggregate by ranks
  ta.list <- list()
  pre.ra <- NA
  for (ra in col.names) {
    ra.col=which(colnames(taxa.assign)==ra)
    
    ### preprocess rank columns
    # BLAST + MEGAN
    id.match <- grep("root|cellular organisms|No hits|Not assigned|unclassified sequences", 
                     taxa.assign[, ra.col], ignore.case = TRUE)
    if (length(id.match) > 0)
      taxa.assign[id.match, ra.col] <- paste("unclassified")
    
    if (length(ta.list) > 0) {
      # replace repeated high rank taxa to unclassified high rank
      # MEGAN unclassified
      id.match <- intersect(!grep("unclassified", taxa.assign[, pre.ra.col], ignore.case = T), 
                            which(taxa.assign[, pre.ra.col]==taxa.assign[, ra.col]))
      # MEGAN environmental samples
      id.match <- c(id.match, grep("environmental samples", taxa.assign[, ra.col], ignore.case = T))
      # RDP unclassified
      id.match <- c(id.match, which(trimAll(taxa.assign[, ra.col])==""))
      if (length(id.match) > 0)
        taxa.assign[id.match, ra.col] <- paste("unclassified", taxa.assign[id.match, pre.ra.col])
      
      # replace unclassified ??? to unclassified high rank
      id.match <- grep("unclassified ", taxa.assign[, pre.ra.col], ignore.case = TRUE)
      if (length(id.match) > 0)
        taxa.assign[id.match, ra.col] <- paste("unclassified", pre.ra)
    }
    
    ### aggregate
    # taxa.assign 1st col is "row.names", "ncol.samples" columns abundence, and length(col.names) columns rank
    data.ta <- taxa.assign[,c(2:(1+ncol.cm), ra.col)]
    taxa.assign.ra <- aggregate(as.formula(paste(". ~", ra)), data=data.ta, FUN=sum)

    # move "unclassified ???" to last
    taxa.assign.ra <- ComMA::mvRowsToLast(taxa.assign.ra, ra, "unclassified")
    # move "unclassified phylum" to last
    if (!is.na(pre.ra)) 
      taxa.assign.ra <- ComMA::mvRowsToLast(taxa.assign.ra, ra, paste0("^unclassified ", pre.ra, "$"))
    # move "unclassified" to last
    taxa.assign.ra <- ComMA::mvRowsToLast(taxa.assign.ra, ra, "^unclassified$")
    
    ta.list[[ra]] <- taxa.assign.ra
    pre.ra.col <- ra.col
    pre.ra <- ra
  }
  return(ta.list)
}

