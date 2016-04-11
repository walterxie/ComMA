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
#' @param min.conf The confidence threshold to drop rows < \emph{min.conf}.
#' @keywords utils
#' @export
#' @examples 
#' tt.sub <- subsetTaxaTable(tt.megan, taxa.group="Proteobacteria", rank="phylum")
#' tt.sub <- subsetTaxaTable(tt.megan, taxa.group="Cnidaria|Brachiopoda|Echinodermata|Porifera", rank="phylum", include=FALSE)
#' 
#' @rdname utilsTaxa
assignTaxa <- function(community.matrix, taxa.table, classifier="MEGAN", min.conf=0.8, 
                       col.names=c("kingdom", "phylum", "class", "order", "family", "genus")) {
  ranks <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  if (!all(col.names %in% ranks)) 
    stop(cat("Invaild column names for ranks !\nUse a subset of", paste(ranks, collapse = ",")))
  
  if (classifier=="RDP") {
    if (! "confidence" %in% colnames(taxa.table))
      stop("Invaild file format from RDP: column 'confidence' is required !")
    
    # drop rows < min.conf
    n.row <- nrow(df.path)
    df.path <- df.path[df.path[,"confidence"] >= min.conf, ] 
    cat("Drop", n.row-nrow(df.path), "rows whose confidence <", min.conf, ".\n")
  }
  
}

