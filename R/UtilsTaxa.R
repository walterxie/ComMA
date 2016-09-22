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
#' @param rank The rank to specify which column name in \code{taxa.table} to search. 
#' @param include Define whether include or exclude given \code{taxa.group}. 
#' Default to TRUE.
#' @param verbose More details. Default to TRUE.
#' @keywords taxonomy
#' @export
#' @examples 
#' tt.sub <- subsetTaxaTable(tt.megan, taxa.group="Proteobacteria", rank="phylum")
#' tt.sub <- subsetTaxaTable(tt.megan, taxa.group="Cnidaria|Brachiopoda|Echinodermata|Porifera", rank="phylum", include=FALSE)
#' 
#' @rdname utilsTaxa
subsetTaxaTable <- function(taxa.table, taxa.group="assigned", rank="kingdom", include=TRUE) {
  # get attr if taxa.table is cm.taxa
  attrs <- attributes(taxa.table)
  ncol.cm <- attrs$ncol.cm
  col.ranks <- attrs$col.ranks
  
  if (include) {
    # include PROTISTS, taxa.group="CHROMISTA|PROTOZOA", rank="kingdom"
    taxa.table <- subset(taxa.table, (grepl(taxa.group, taxa.table[,rank], ignore.case = T))) 
  } else { 
    # exclude some phyla, taxa.group="Cnidaria|Brachiopoda|Echinodermata|Porifera", rank="phylum"
    taxa.table <- subset(taxa.table, !grepl(taxa.group, taxa.table[,rank], ignore.case = T)) 
  }
  
  # add attr if they are not NULL
  if (! is.null(ncol.cm))
    attr(taxa.table, "ncol.cm") <- ncol.cm
  if (! is.null(col.ranks)) 
    attr(taxa.table, "col.ranks") <- col.ranks
  
  return(taxa.table)
}

#' @details \code{subsetCM} returns a subset community matrix 
#' regarding \code{taxa.group}. 
#' Set either \code{taxa.group} or \code{rank} to NA, as default, 
#' to use the whole \code{taxa.table}.
#' 
#' @keywords taxonomy
#' @export
#' @examples 
#' sub.cm <- subsetCM(cm, tt, taxa.group="BACTERIA", rank="kingdom")
#' 
#' @rdname utilsTaxa 
subsetCM <- function(community.matrix, taxa.table, taxa.group=NA, rank=NA, 
                     col.ranks=c("kingdom", "phylum", "class", "order", "family")) {
  if (is.na(taxa.group) || is.na(rank))
    tt.sub <- taxa.table
  else
    tt.sub <- ComMA::subsetTaxaTable(taxa.table, taxa.group=taxa.group, rank=rank)	
  cm.taxa <- ComMA::mergeCMTaxa(community.matrix, tt.sub, col.ranks=col.ranks, has.total=0)
  cm.taxa <- cm.taxa[,colnames(community.matrix)]
  return(cm.taxa)
}

#' @details \code{prepTaxonomy} replace repeated high rank taxa to 
#' unclassified high rank in MEGAN result, 
#' or replace the blank value to unclassified in RDP result,
#' in order to make taxonomy table \code{taxa.table} (can be \code{cm.taxa}) 
#' to make names look nice.
#' \code{col.ranks} vector have to be rank column names in taxa.table.
#' 
#' @param pattern The pattern for \code{\link{gsub}} "perl = TRUE". 
#' Default to "(\\s\\[|\\()(\\=|\\.|\\,|\\s|\\w|\\?)*(\\]|\\))".
#' @param txt.unclassified The key word to represent unclassified taxonomy.
#' @keywords taxonomy
#' @export
#' @examples 
#' tt <- prepTaxonomy(taxa.table, col.ranks=c("kingdom", "phylum", "class"))
#' 
#' @rdname utilsTaxa 
prepTaxonomy <- function(taxa.table, col.ranks=c("kingdom", "phylum", "class", "order", "family"),
                         txt.unclassified="unclassified",
                         pattern="(\\s\\[|\\()(\\=|\\.|\\,|\\s|\\w|\\?)*(\\]|\\))") {
  parent.rank <- NA
  for (rank in col.ranks) {
    if (! rank %in% colnames(taxa.table))
      stop("Invalid ", rank, " not existing in colnames(taxa.table) !\n", 
           paste(colnames(taxa.table), collapse = ","))
    
    taxa.table[, rank] <- gsub("root|cellular organisms|No hits|Not assigned|unclassified sequences", 
                               txt.unclassified, taxa.table[, rank], ignore.case = TRUE)
    # Remove assorted quirks in taxonomy! 
    taxa.table[, rank] <- gsub(pattern, "", taxa.table[, rank], perl = TRUE)
    
    # MEGAN unclassified
    if (!is.na(parent.rank)) {
      ta.tmp <- subset(taxa.table, !grepl("unclassified", taxa.table[, parent.rank], ignore.case = T))
      id.match <- which(tolower(ta.tmp[, parent.rank])==tolower(ta.tmp[, rank]))
      
      # MEGAN environmental samples => unclassified ???
      id.match <- c(id.match, grep("environmental samples", taxa.table[, rank], ignore.case = T))
      
      # RDP unclassified
      id.match <- c(id.match, which(trimSpace(taxa.table[, rank])==""))
      if (length(id.match) > 0)
        taxa.table[id.match, rank] <- paste("unclassified", taxa.table[id.match, parent.rank])
      
      # replace unclassified ??? to unclassified high rank
      id.match <- grep("unclassified ", taxa.table[, parent.rank], ignore.case = TRUE)
      if (length(id.match) > 0)
        taxa.table[id.match, rank] <- paste("unclassified", parent.rank)
    } 
    parent.rank <- rank
  }
  return(taxa.table)
}

#' @details 
#' \code{mergeCMTaxa} creates a data frame \code{cm.taxa} combined community matrix with 
#' taxonomic classification table. The 1st column is "row.names" that are OTUs/individuals, 
#' the next "ncol.cm" columns are abundence that can be sample-based or total, 
#' and the last "length(col.ranks)" columns are the ranks. 
#' 
#' All sequences either classified as 
#' "root|cellular organisms|No hits|Not assigned|unclassified sequences"
#' from BLAST + MEGAN, or confidence < \emph{min.conf} threshold from RDP, 
#' are changed to "unclassified", which will be moved to the last row.  
#' 
#' @param community.matrix Community matrix (OTU table), where rows are 
#' OTUs or individual species and columns are sites or samples. See \code{\link{ComMA}}.
#' @param classifier The classifier is used to generate \code{taxa.table}. 
#' Value is MEGAN or RDP. Default to MEGAN.
#' @param min.conf The confidence threshold to drop rows < \emph{min.conf}.
#' @param has.total If 0, then only return abundance by samples (columns) of community matrix. 
#' If 1, then only return total abundance. If 2, then return abundance by samples (columns) and total. 
#' Default to 1.
#' @param preprocess If TRUE, as default, replace 
#' "root|cellular organisms|No hits|Not assigned|unclassified sequences" from MEGAN result, 
#' or mark OTUs as 'unclassified' in RDP result whose confidence < \code{min.conf} threshold. 
#' @param col.ranks A vector or string of column name(s) of taxonomic ranks in the taxa table, 
#' which will determine the aggregated abundence matrix. They have to be full set or subset of 
#' c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"). 
#' Default to c("kingdom", "phylum", "class", "order", "family").
#' @param sort Sort the taxonomy rank by rank. Default to TRUE.
#' @param mv.row.names Default to TRUE to move the column 'Row.names' 
#' created by \code{\link{merge}} into data frame row.names, 
#' in order to keep the 1st column same as community matrix. 
#' Suggest not to change it.
#' @return 
#' \code{ncol.cm} and \code{col.ranks} are attributes of \code{cm.taxa} generated by \code{mergeCMTaxa}.
#' 
#' \code{ncol.cm} indicates how many column(s) is/are abundence in \code{cm.taxa}.
#'
#' \code{col.ranks} records what ranks column(s) is/are in \code{cm.taxa}, 
#' which is also the input of \code{mergeCMTaxa}. 
#' 
#' @keywords taxonomy
#' @export
#' @rdname utilsTaxa
mergeCMTaxa <- function(community.matrix, taxa.table, classifier=c("MEGAN","RDP"), min.conf=0.8, 
                        has.total=1, sort=TRUE, preprocess=TRUE, verbose=TRUE, 
                        mv.row.names=T, 
                        col.ranks=c("kingdom", "phylum", "class", "order", "family")) {
  classifier <- match.arg(classifier)
  if (verbose)
    cat("Parse classification from", classifier, "classifier.\n")
  ranks <- getRanks()
  if (length(col.ranks) < 1 || !all(col.ranks %in% ranks)) 
    stop("Invaild column names for ranks !\nUse one element or a subset of ", 
         paste(ranks, collapse = ","))
  if (!all(col.ranks %in% colnames(taxa.table))) 
    stop("Column names in taxa.table do not have ", paste(col.ranks, collapse = ","))
  if (! "kingdom" %in% colnames(taxa.table))
    stop("Column names in taxa.table must have 'kingdom' column !")
  if ("confidence" %in% colnames(taxa.table) && classifier != "RDP")
    stop("Find 'confidence' column, please define classifier = 'RDP' !")
  
  ###### merge community.matrix, taxa.table 
  cm <- data.frame(row.names = rownames(community.matrix))
  if (has.total!=1) 
    cm <- community.matrix
  if (has.total > 0) 
    cm[,"total"] <- rowSums(community.matrix) 
  ncol.cm <- ncol(cm)
  
  if (classifier == "RDP") {
    if (! "confidence" %in% colnames(taxa.table))
      stop("Invaild file format from RDP: column 'confidence' is required !")
    
    ### grep all "Unclassified" in kingdom from greengenes, and overwrite to "unclassified"
    id.match <- grep("^Unclassified$", taxa.table[, "kingdom"], ignore.case = TRUE)
    if (length(id.match) > 0) {
      taxa.table[id.match, col.ranks] <- "unclassified"
      if (verbose)
        cat("Find", length(id.match), "rows marked as 'unclassified' from taxa table.\n")
    }
    
    ### set rows < min.conf as 'unclassified'
    n.row.un <- nrow(taxa.table[taxa.table[,"confidence"] < min.conf, ])
    taxa.table[taxa.table[,"confidence"] < min.conf, col.ranks] <- "unclassified"
    
    if (verbose)
      cat("Set", n.row.un, "rows as 'unclassified' from the total of", nrow(taxa.table), 
          "in RDP taxa table, whose confidence <", min.conf, ".\n")
  } # else BLAST + MEGAN

  if (preprocess)
    taxa.table <- prepTaxonomy(taxa.table, col.ranks=col.ranks)
  
  # cm.taxa 1st col is "row.names", "ncol.cm" columns abundence, and length(col.ranks) columns rank
  cm.taxa <- merge(cm, taxa.table, by = "row.names")
  # move 'Row.names' into row.names to keep 1st column same as community matrix.
  if (mv.row.names) {
    rownames(cm.taxa) <- cm.taxa[,"Row.names"]
    cm.taxa <- cm.taxa[,-1]
  }
  
  if (sort) {
    # order rank by rank
    ord.cmd = parse(text = paste0('cm.taxa[order(cm.taxa[,"', paste(col.ranks, collapse = '"], cm.taxa[,"'), '"]),]')) 
    cm.taxa <- eval(ord.cmd)
  }
  
  if (verbose)
    cat("Merge", nrow(cm), "rows in community matrix with", nrow(taxa.table), "rows in taxa table, get", 
        nrow(cm.taxa), "classified OTUs.\n")
  
  attr(cm.taxa, "ncol.cm") <- ncol.cm
  attr(cm.taxa, "col.ranks") <- col.ranks
  return(cm.taxa)
}

#' @details 
#' \code{assignTaxaByRank} provides a list of taxonomic assignments with abundance 
#' from community matrix at different rank levels, where rownames are taxonomy 
#' at that rank, and columns are the sample names (may include total). 
#' The function is iterated through \code{col.ranks}, and \code{\link{aggregate}}s 
#' abundance into taxonomy based on the rank in \code{col.ranks}. 
#' 
#' @param cm.taxa The data frame combined community matrix with 
#' taxonomic classifications generated by \code{mergeCMTaxa}.
#' The row.names are OTUs, 1st column is the start of community matrix, 
#' \code{ncol.cm} column is the end of abundence, 
#' and \code{length(col.ranks)} columns taxonomy at different ranks. 
#' It should have attributes \code{ncol.cm} and \code{col.ranks}.
#' 
#' Note: From 1 to \code{ncol.cm} columns, the last column may be 'total' 
#' that is rowSums(cm) determined by \code{has.total} in \code{mergeCMTaxa}.
#' 
#' @param unclassified An interger to instruct how to deal with 
#' "unclassified" taxa. Default to 0, which keeps all "unclassified"
#' but moves them to the last rows. 
#' If 1, then remove the row whose taxon name is exact "unclassified". 
#' See the detail. 
#' If 2, then remove the row whose taxon name is exact "unclassified", 
#' but also merge all the rest "unclassified ???" to "unclassified rank",
#' such as "unclassified family".
#' If 3, then remove every rows containing "unclassified".
#' if 4, then do nothing.
#' @param aggre.FUN A function for \code{FUN} in \code{\link{aggregate}}. 
#' Default to \code{sum} to provide the reads abundance. 
#' Make \code{aggre.FUN=function(x) sum(x>0)} provide the OTU abundance.
#' @keywords taxonomy
#' @export
#' @examples 
#' cm.taxa <- mergeCMTaxa(community.matrix, tt.megan) 
#' ta.megan <- assignTaxaByRank(cm.taxa)
#' 
#' cm.taxa <- mergeCMTaxa(community.matrix, tt.rdp, classifier="RDP", has.total=0)
#' ta.rdp <- assignTaxaByRank(cm.taxa, unclassified=2)
#' colSums(ta.rdp[["phylum"]])
#' 
#' @rdname utilsTaxa
assignTaxaByRank <- function(cm.taxa, unclassified=0, aggre.FUN=sum) {
  attr.cm.ta <- attributes(cm.taxa)
  ncol.cm <- attr.cm.ta$ncol.cm
  col.ranks <- attr.cm.ta$col.ranks
  
  if (is.null(ncol.cm) || is.null(col.ranks)) 
    stop("Input cm.taxa should have attributes ncol.cm and col.ranks !")
  
  ### preprocess rank columns
  if (unclassified < 4)
    cm.taxa <- prepTaxonomy(cm.taxa, col.ranks=col.ranks)
  
  ###### aggregate by ranks
  ta.list <- list()
  pre.ra <- NA
  for (ra in col.ranks) {
    ra.col=which(colnames(cm.taxa)==ra)
    
    ### aggregate
    # cm.taxa 1st col is "row.names", "ncol.samples" columns abundence, and length(col.ranks) columns rank
    data.ta <- cm.taxa[,c(1:ncol.cm, ra.col)]
    # e.g. family plot1 plot2 ...
    taxa.assign <- aggregate(as.formula(paste(". ~", ra)), data=data.ta, FUN=aggre.FUN)

    # deal with "unclassified"
    if (unclassified==0) {
      # move "unclassified ???" to last
      taxa.assign <- ComMA::mvRowsToLast(taxa.assign, ra, "unclassified")
      # move "unclassified phylum" to last
      if (!is.na(pre.ra)) 
        taxa.assign <- ComMA::mvRowsToLast(taxa.assign, ra, paste0("^unclassified ", pre.ra, "$"))
      # move "unclassified" to last
      taxa.assign <- ComMA::mvRowsToLast(taxa.assign, ra, "^unclassified$")
    } else if (unclassified==1) {
      taxa.assign <- subset(taxa.assign, !grepl("^unclassified$", taxa.assign[, ra], ignore.case = T))
    } else if (unclassified==2) {
      ta.tmp <- subset(taxa.assign, !grepl("unclassified", taxa.assign[, ra], ignore.case = T))
      ta.un.tmp <- subset(taxa.assign, grepl("unclassified", taxa.assign[, ra], ignore.case = T))
      ta.un.tmp <- subset(ta.un.tmp, !grepl("^unclassified$", ta.un.tmp[, ra], ignore.case = T))
      if (nrow(ta.un.tmp) > 0) {
        # rank + colnames(cm)
        if (ncol(ta.un.tmp) != ncol.cm + 1)
          stop("Invalid index, unclassified =", unclassified, ", where ncol", 
               ncol(ta.un.tmp), "should ==", ncol.cm + 1, "!")
        
        if (ncol(ta.un.tmp) > 2)
          ra.un <- c(paste("unclassified", ra), colSums(ta.un.tmp[,-1]))
        else
          ra.un <- c(paste("unclassified", ra), sum(ta.un.tmp[,-1]))
        ta.tmp <- rbind(ta.tmp, ra.un)
      }
      taxa.assign <- ta.tmp
    } else if (unclassified==3) {
      taxa.assign <- subset(taxa.assign, !grepl("unclassified", taxa.assign[, ra], ignore.case = T))
    }
    taxa.assign[,2:ncol(taxa.assign)] <- sapply(taxa.assign[,2:ncol(taxa.assign)], as.numeric)
    
    # move taxa to rownames
    rownames(taxa.assign) <- taxa.assign[,1]
    # add rank attr
    attr(taxa.assign, "rank") <- ra
    # make sure 1-column data frame not converted to vector, see ?"[" 
    taxa.assign <- taxa.assign[,-1, drop=FALSE]
    
    # rm 0 total row
    ta.list[[ra]] <- taxa.assign
    pre.ra.col <- ra.col
    pre.ra <- ra
  }
  return(ta.list)
}

#' @details 
#' \code{groupsTaxaMembers} groups the members (rows, also OTUs) from 
#' \code{cm.taxa} for each taxa in \code{taxa.assign} at the \code{rank}, 
#' and returns a list of members (OTUs) grouped by taxonomy. 
#' Default to drop all unclassified members (OTUs). 
#' 
#' It is impossible to trace back members after \code{assignTaxaByRank},
#' so that this function only has one option except the default, 
#' which assign the rest of members (OTUs) not picked up from other taxa  
#' into "unclassified". The result relies on using the identical \code{cm.taxa} 
#' in both \code{assignTaxaByRank} and \code{groupsTaxaMembers}.
#' 
#' @param taxa.assign The data frame of taxonomic assignments with abundance
#' at the \code{rank}, where rownames are taxonomy at that rank, 
#' and columns are the sample names (may include total). It can be 
#' one element of the list generated by \code{assignTaxaByRank}. 
#' See the detail.
#' @param regex1,regex2 Use for \code{\link{gsub}(regex1, regex2, row.names)} 
#' to remove or replace annotation from original labels. 
#' Default to \code{regex1="(\\|[0-9]+)", regex2=""}, 
#' which removes size annotation seperated by "|".
#' @param rm.unclassified Drop all unclassified rows (OTUs). Default to TRUE.
#' @keywords taxonomy
#' @export
#' @examples 
#' taxa.members <- groupsTaxaMembers(ta.rdp[["phylum"]], tt.rdp)
#' taxa.members <- groupsTaxaMembers(ta.rdp[["family"]], tt.rdp, rank="family")
#' 
#' @rdname utilsTaxa
groupsTaxaMembers <- function(taxa.assign, cm.taxa, rank="phylum", rm.unclassified=TRUE,
                              regex1="(\\|[0-9]+)", regex2="", ignore.case=TRUE, 
                              verbose=TRUE) {
  if (is.null(rownames(taxa.assign)) || is.null(rownames(cm.taxa)))
    stop("Invalid taxa.assign or cm.taxa, which needs taxa or OTUs in its rownames !")
  ra <- attributes(taxa.assign)$rank
  if (is.null(ra))
    warning("Cannot find rank attribute from taxa.assign, no validation is proceeded !")
  else if (ra != rank)
    stop("Inconsistent rank level for taxa.assign : ", ra, " != ", rank, " !")
  
  id.unclassified <- grep("unclassified", rownames(taxa.assign), ignore.case = T)
  if (length(id.unclassified) > 0) {
    if (rm.unclassified)
      warning("Drop all unclassified rows (OTUs) belonging to : ", 
              paste(rownames(taxa.assign)[id.unclassified], collapse = ","), " !\n")
    
    taxa.assign <- taxa.assign[-id.unclassified,]
  }
  
  # remove/replace annotation
  if (! is.null(regex1)) 
    rownames(cm.taxa) <- gsub(regex1, regex2, rownames(cm.taxa), ignore.case = ignore.case)
  
  taxa.members <- sapply(rownames(taxa.assign), 
                         function(taxa) rownames(cm.taxa)[which(cm.taxa[,rank] == taxa)])
  
  # add the rest of OTUs into unclassified
  if (!rm.unclassified && length(id.unclassified) > 0) {
    uncl.members <- setdiff(rownames(cm.taxa), unlist(taxa.members)) 
    taxa.members[["unclassified"]] <- uncl.members
    if (verbose)
      cat("Add", length(uncl.members), "OTUs into unclassified.\n")
  }
  
  if (!all(lapply(taxa.members,length)>0)) {
    warning("Remove taxa having no OTU : ", 
            paste(names(taxa.members)[lapply(taxa.members,length)<=0], collapse = ","), " !\n")
    # Remove empty elements from list
    taxa.members <- taxa.members[lapply(taxa.members,length)>0]
  }
  
  return(taxa.members)
}

# taxonomic ranks for column names
getRanks <- function() {
  ranks <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
}

