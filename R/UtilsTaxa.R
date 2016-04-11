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
#' \code{assignTaxa} provides taxonomic assignment with abudence from community matrix.
#' 
#' @param community.matrix Community matrix (OTU table), where rows are 
#' OTUs or individual species and columns are sites or samples. See \code{\link{ComMA}}.
#' @param min.conf The confidence threshold to drop rows < \emph{min.conf}.
#' @param total If 0, then only return abudence by samples (columns) of community matrix. 
#' If 1, then only return toal abudence. If 2, then return abudence by samples (columns) and total. 
#' Default to 1.
#' @keywords utils
#' @export
#' @examples 
#' tt.sub <- subsetTaxaTable(tt.megan, taxa.group="Proteobacteria", rank="phylum")
#' tt.sub <- subsetTaxaTable(tt.megan, taxa.group="Cnidaria|Brachiopoda|Echinodermata|Porifera", rank="phylum", include=FALSE)
#' 
#' @rdname utilsTaxa
assignTaxa <- function(community.matrix, taxa.table, classifier="MEGAN", min.conf=0.8, total=1,
                       col.names=c("kingdom", "phylum", "class", "order", "family", "genus")) {
  ranks <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  ###### preprocess rank columns
  
  
  ###### merge community.matrix, taxa.table 
  if (!all(col.names %in% ranks)) 
    stop(cat("Invaild column names for ranks !\nUse a subset of", paste(ranks, collapse = ",")))

  cm <- data.frame(row.names = rownames(community.matrix))
  if (total!=1) 
    cm <- community.matrix
  if (total > 0) 
    cm$Total <- rowSums(community.matrix) 
  
  if (classifier=="RDP") {
    if (! "confidence" %in% colnames(taxa.table))
      stop("Invaild file format from RDP: column 'confidence' is required !")
    
    # drop rows < min.conf
    n.row <- nrow(taxa.table)
    taxa.table <- taxa.table[taxa.table[,"confidence"] >= min.conf, ] 
    cat("Drop", n.row-nrow(taxa.table), "rows from the total of", n.row, 
        "in RDP taxa table, whose confidence <", min.conf, ".\n")
  }
  
  taxa.assign <- merge(cm, taxa.table, by = "row.names")
  
  cat("Merge", nrow(cm), "rows in community matrix with", nrow(taxa.table), "rows in taxa table, get", 
      nrow(taxa.assign), "classifications.\n")
  
  ###### aggregate by ranks
  
}

