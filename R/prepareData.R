
#' @name getData
#' @title Get example data for \pkg{ComMA} package
#'
#' @description Data \strong{communityMatrix} \strong{t.communityMatrix} is 
#' a transposed matrix from community matrix we defined here.
#' \strong{phyloTree} is a rooted tree of phylo object, 
#' which can get from \pkg{ape} \code{\link{read.tree}}.
#' 
# Add postfix for figure name or table label
#' @param ... A list of names to be concatenated
#' @param sep Default to dash "-"
#' @return \code{postfix} returns a name concatenated by substrings and separated by \code{sep}.
#' @export
#' @examples 
#' postfix("16S", "assigned")
#'[1] "16S-assigned-byplot-min2"
#'
#' @rdname getData
postfix <- function(..., isPlot=TRUE, taxa.group="all", minAbund=1, minRich=0, sep="-") {
  full.name <- paste(list(...), collapse = sep)
  
  if (taxa.group!="all") 
    full.name <- paste(full.name, "taxa.group", sep = sep) 
  if (isPlot) 
    full.name <- paste(full.name, "byplot", sep = sep) 
  if (minAbund > 0) 
    full.name <- paste(full.name, paste0("min", minAbund+1), sep = sep)
  if (minRich > 0) 
    full.name <- paste(full.name, minRich+1, sep = sep)
  
  return(full.name)
}

#' Extract plot names from full names (plot + subplot) separated by \code{sep}.
#' 
#' @param full.name The full name has plot and subplot together, but separated by \code{sep}.
#' @return \code{getPlot} returns plot names.
#' @export
#' @examples 
#' getPlot(c("Plot1-A", "CM30c39-L"))
#' [1] "Plot1" "CM30c39"
#' 
#' @rdname getData
getPlot <- function(full.name, sep="-") {
  sapply(strsplit(as.character(full.name), sep), "[[", 1)
}

######## load community matrix #######
# expId = 1:6
# isPlot determines to use which matrix file, by subplot or plot 
# min2 = rmSingleton, whether remove all singletons
# may contain empty rows cols
#' @param matrix.name The string to locate the matrix from its file name.
#' @param isPlot Boolean value to determine the matrix file sampled by subplot or plot
#' @param minAbund The minimum abundance threshold to determine which 
#' row/column to be removed. For exampe, if minAbund=1, then remove 
#' all singletons appeared in only one sample. Default to 1 (singletons).
#' But \code{postfix} is only used for naming, no data processed.
#' @param verbose More details. Default to TRUE.
#' @return \code{getCommunityMatrix} returns a community matrix, 
#' where rows are OTUs or individual species and columns are sites or samples. 
#' It may contain empty rows or columns. \code{prepCommunityMatrix} can be used 
#' to remove them.
#' @export
#' @examples 
#' communityMatrix <- getCommunityMatrix("16S", isPlot=TRUE, minAbund=0)
#' 
#' @rdname getData
getCommunityMatrix <- function(matrix.name, isPlot, minAbund=1, verbose=TRUE) {
  if (isPlot) {
    inputCM <- file.path("data", paste(matrix.name, "by_plot.txt", sep="_"))
  } else {
    # e.g. data/16S.txt
    inputCM <- file.path("data", paste(matrix.name, ".txt", sep=""))
  }
  
  communityMatrix <- readFile(inputCM)
  if(verbose) 
    cat("\nUpload community matrix : ", ncol(communityMatrix), "columns,", nrow(communityMatrix), "rows, from", inputCM, "\n") 
  
  rmMinAbundance(communityMatrix, minAbund)
  
  return(communityMatrix)
}

#' remove empty rows or columns after \code{getCommunityMatrix}.
#' @return \code{prepCommunityMatrix} returns a community matrix no empty rows or columns.
#' @export
#' @examples 
#' communityMatrix <- prepCommunityMatrix(communityMatrix)
#' 
#' @rdname getData
prepCommunityMatrix <- function(communityMatrix) {
  # filter column first to avoid empty rows after columns remvoed if vectorThr>0
  if(any(colSums(communityMatrix)== 0)) {
    communityMatrix <- rmVectorFromCM(communityMatrix, vectorThr=0, MARGIN=2)
    #stop("Invalid input: community matrix has empty samples !")
  }
  if(any(rowSums(communityMatrix)== 0)) {
    communityMatrix <- rmVectorFromCM(communityMatrix, vectorThr=0, MARGIN=1)
    #stop("Invalid input: community matrix has empty OTUs !")
  }
  
  return(communityMatrix)
}

#' Get a transposed matrix of community matrix, given more filters, such as \code{taxa.group}, \code{minRich}.
#' 
#' @param minRich The minimum richness to keep matrix. For example, 
#' drop the matrix (return NULL) if OTUs less than this threshold. Default to 0.
#' @return \code{getCommunityMatrixT} returns a transposed matrix of community matrix, 
#' where columns are OTUs or individual species and rows are sites or samples. 
#' It is also the abundances argument in \pkg{vegetarian} \code{\link{d}}.
#' return \emph{NULL}, if no OTUs or OTUs less than \code{minRich} threshold.
#' @export
#' @examples 
#' t.communityMatrix <- getCommunityMatrixT("16S", isPlot=TRUE, minAbund=1, taxa.group="BACTERIA")
#' 
#' @rdname getData
getCommunityMatrixT <- function(matrix.name, isPlot, taxa.group="all", minAbund=1, minRich=0, verbose=TRUE) {
  communityMatrix <- getCommunityMatrix(matrix.name, isPlot, minAbund)
  
  if (taxa.group != "all") {
    ##### load data #####
    taxaPaths <- getTaxaPaths(matrix.name, taxa.group)
    
    if (nrow(taxaPaths) < minRich) {
      cat("Warning: return NULL, because", nrow(taxaPaths), "row(s) classified as", taxa.group, "<", minRich, "threshold.\n")
      return(NULL)
    } else {
      # merge needs at least 2 cols 
      taxaAssgReads <- merge(communityMatrix, taxaPaths, by = "row.names")
      
      if (nrow(taxaAssgReads) < minRich) {
        cat("Warning: return NULL, because", nrow(taxaAssgReads), "row(s) in cm match", taxa.group, "<", minRich, "threshold.\n")
        return(NULL)
      }
      
      # move 1st col Row.names to row.names
      rownames(taxaAssgReads) <- taxaAssgReads[,"Row.names"]
      taxaAssgReads <- taxaAssgReads[,-1]
      # get CM
      taxaAssgReads <- taxaAssgReads[,1:ncol(communityMatrix)]
      
      cat("Merging", nrow(taxaAssgReads), "matched OTUs from", nrow(communityMatrix), "OTUs in matrix to", 
          nrow(taxaPaths), "taxa classification, taxa.group =", taxa.group, ", final ncol =", ncol(taxaAssgReads), ".\n")
      
      communityMatrix <- data.matrix(taxaAssgReads)
    }
  }
  
  communityMatrix <- prepCommunityMatrix(communityMatrix)
  
  communityMatrixT <- transposeCM(communityMatrix)
  
  return(communityMatrixT)
}
# rownames(communityMatrix) <- gsub("-(.*)|", "\\1", rownames(communityMatrix))

###### taxa assignment by reads #####
# "ARCHAEA", "BACTERIA", "CHROMISTA", "PROTOZOA", "FUNGI", "PLANTAE", "ANIMALIA", "EUKARYOTA", "PROKARYOTA", "PROTISTS"
# PROKARYOTA: all prokaryotes (Bacteria + Archaea)
# EUKARYOTA: all eukaryotes
# PROTISTS: "CHROMISTA|PROTOZOA" all micro-eukaryotes
#' @param taxa.group The taxonomic group, the values can be 'all', 'assigned', or 
#' Group 'all' includes everything.
#' Group 'assigned' removes all uncertain classifications including 
#' 'root', 'cellular organisms', 'No hits', 'Not assigned'. 
#' Alternatively, any high-ranking taxonomy in your taxonomy file 
#' can be used as a group, such as 'BACTERIA', 'Proteobacteria', etc.
#' @return \code{getCommunityMatrixT} returns a transposed matrix of community matrix, 
#' where columns are OTUs or individual species and rows are sites or samples. 
#' It is also the abundances argument in \pkg{vegetarian} \code{\link{d}}.
#' return \emph{NULL}, if no OTUs or OTUs less than \code{minRich} threshold.
#' @export
#' @examples 
#' taxaPaths <- getTaxaPaths(matrix.name="18S", taxa.group="assigned")
#' 
#' @rdname getData
getTaxaPaths <- function(matrix.name, taxa.group="all", rank="kingdom", verbose=TRUE) {
  inputTaxa <- file.path("data", "taxonomy_tables", paste(matrix.name, "taxonomy_table.txt", sep="_"))
  taxaPaths <- readTaxaFile(inputTaxa)	
  taxaPaths <- taxaPaths[order(rownames(taxaPaths)),]
  # make lower case to match ranks
  colnames(taxaPaths) <- tolower(colnames(taxaPaths))
  
  nTaxa=nrow(taxaPaths)
  ##### keep OTU rows contain given taxa belongTo ##### 
  if (taxa.group != "all") {
    # Exclude unassigned etc
    taxaPaths <- subset(taxaPaths, !(grepl("root|cellular organisms|No hits|Not assigned", taxaPaths[,"kingdom"])))  
    # (Only retain prokaryotes for 16S, eukaryotes for the other amplicons)  
    if (toupper(matrix.name)=="16S") {
      taxaPaths <- subset(taxaPaths, (grepl("BACTERIA|ARCHAEA", taxaPaths[,"kingdom"])))  
    } else {
      taxaPaths <- subset(taxaPaths, !(grepl("BACTERIA|ARCHAEA", taxaPaths[,"kingdom"])))
    }
    # Exclude probably bogus taxa
    taxaPaths <- subset(taxaPaths, !(grepl("Cnidaria|Brachiopoda|Echinodermata|Porifera", taxaPaths[,"phylum"])|
                                       grepl("Bivalvia|Teleostei|Elasmobranchii|Polyplacophora", taxaPaths[,"class"])|
grepl("Nudibranchia|Crocodylia|Serpentes|Testudines|Carnivora|Gymnophiona|Lagomorpha|Rodentia|Serpentes|Scorpiones", taxaPaths[,"order"])))
    
    if (taxa.group != "assigned") {
      if (toupper(taxa.group) == "PROKARYOTA" || toupper(taxa.group) == "EUKARYOTA") {
        taxaPaths <- subset(taxaPaths, grepl(taxa.group, taxaPaths[,"superkingdom"])) 
      } else if (toupper(taxa.group) == "PROTISTS") {
        taxaPaths <- subset(taxaPaths, grepl("CHROMISTA|PROTOZOA", taxaPaths[,"kingdom"])) 
      } else {
        taxaPaths <- subset(taxaPaths, grepl(taxa.group, taxaPaths[,rank])) 
      }
    }
  }
  
  if (nrow(taxaPaths) < 1)
    cat("Warning: cannot find", taxa.group, "from taxa path file", inputTaxa, "!")
  
  if(verbose && nrow(taxaPaths) < nTaxa) 
    cat("\nSelect", nrow(taxaPaths), "taxa classification, taxa.group =", taxa.group, ".\n") 
  
  return(taxaPaths)
}


# rankLevel: the taxa level in each bar
# groupLevel: used to assign colour for each group, and must higher than rankLevel
# taxa.group: keep OTU rows contain given taxa group, if "all", keep all
# return taxaAssgReads = CM + rankLevel + groupLevel
getTaxaAssgReads <- function(matrix.name, isPlot, minAbund=1, rankLevel, groupLevel, taxa.group="all") {
  cat("Create taxonomy assignment for", matrix.name, ".\n")
  
  ##### load data #####
  communityMatrix <- getCommunityMatrix(matrix.name, isPlot, minAbund)
  communityMatrix <- prepCommunityMatrix(communityMatrix)
  
  communityMatrix <- communityMatrix[order(rownames(communityMatrix)),]
  communityMatrix <- communityMatrix[,order(colnames(communityMatrix))]
  
  taxaPaths <- getTaxaPaths(matrix.name, taxa.group)
  
  ###### taxa assignment by reads #####
  if ( ! tolower(rankLevel) %in% tolower(colnames(taxaPaths)) ) 
    stop( paste("Column name", rankLevel, "not exist in taxa path file for ", matrix.name) )
  if (! tolower(groupLevel) %in% tolower(colnames(taxaPaths)) ) 
    stop( paste("Column name", groupLevel, "not exist in taxa path file for ", matrix.name) )
  
  colRankLevel <- which(tolower(colnames(taxaPaths))==tolower(rankLevel))
  colGroupLevel <- which(tolower(colnames(taxaPaths))==tolower(groupLevel))
  
  taxaAssgReads <- merge(communityMatrix, taxaPaths[,c(colRankLevel, colGroupLevel)], by = "row.names")
  
  taxaAssgReads[,rankLevel] <- gsub("(\\s\\[=.*\\])", "", taxaAssgReads[,rankLevel])
  taxaAssgReads[,groupLevel] <- gsub("(\\s\\[=.*\\])", "", taxaAssgReads[,groupLevel])
  
  cat("Merging:", nrow(taxaAssgReads), "OTUs are matched from", nrow(communityMatrix), "OTUs in matrix to", 
      nrow(taxaPaths), "taxa classification, taxa.group =", taxa.group, ".\n")
  
  return(taxaAssgReads)
}

getTaxaRef <- function() {
  tax_ref <- read.table(file.path("data", "New_taxonomy_from_PLOSONE_2015_fixed.txt"), 
                        header = TRUE, sep = "\t", quote = "", comment.char = "")
  # make lower case to match ranks
  colnames(tax_ref) <- tolower(colnames(tax_ref))
  
  # Remove quirks/questions in taxa ([= ...])
  tax_ref <- apply(tax_ref, 2, function(col) gsub("(\\s\\[=.*\\])", "", col))
}

###### Trees #####
getPhyloTree <- function(fNameStem, verbose=TRUE) {
  inputT <- file.path("data", "trees", paste(fNameStem, "tre", sep = "."))
  if (file.exists(inputT)) {
    cat("Load tree from", inputT, "\n") 
    tree <- read.tree(inputT)
    if(verbose) print(tree)
  } else {
    tree <- NULL
    cat("Warning: cannot find tree file:", inputT, "\n") 
  }
  tree
}

###### table to plot Phylo Rarefaction ##### 
getPhyloRareTable <- function(expId, isPlot, min2, taxa.group="assigned", verbose=TRUE) {
  n <- length(matrixNames) 
  # hard code for Vegetation that only has plot and always keep singletons
  if (expId == n) {
    mid.name <- postfix("all", TRUE, FALSE, sep="-")
  } else {
    mid.name <- postfix(taxa.group, isPlot, min2, sep="-") 
  }
  
  inputT <- file.path("data", "pdrf", paste(matrixNames[expId], mid.name, "phylorare", "table.csv", sep="-"))
  if (file.exists(inputT)) {
    phylo.rare.df <- read.csv(file=inputT, head=TRUE, sep=",", row.names=1, check.names=FALSE)
    if(verbose) 
      cat("\nUpload phylo rarefaction table from", inputT, "\n") 
  } else {
    phylo.rare.df <- NULL
    cat("Warning: cannot find phylo rarefaction table", inputT, "\n") 
  }
  phylo.rare.df
}

###### table to plot Rarefaction ##### 
getRarefactionTableTaxa <- function(expId, isPlot, min2, taxa.group, div="alpha1", verbose=TRUE) {
  pathFileStem <- file.path("data", "rf", paste(matrixNames[expId], 
                    postfix(taxa.group, isPlot, rmSingleton, sep="-"), sep = "-"))
  inputT <- paste(pathFileStem, "rare", div, "table.csv", sep = "-")
  if (file.exists(inputT)) {
    rare.df <- read.csv(file=inputT, head=TRUE, sep=",", row.names=1, check.names=FALSE)
    if(verbose) 
      cat("\nUpload rarefaction table per sample from", inputT, "\n") 
  } else {
    rare.df <- NULL
    cat("Warning: cannot find rarefaction table per sample", inputT, "\n") 
  }
  rare.df
}

getRarefactionTable <- function(expId, isPlot, min2, verbose=TRUE) {
  n <- length(matrixNames) 
  matrixName <- matrixNames[expId]
  # hard code for Vegetation that only has plot and always keep singletons
  if (expId == n) {
    matrixName <- postfix(matrixName, TRUE, FALSE, sep="-")
  } else {
    matrixName <- postfix(matrixName, isPlot, min2, sep="-") 
  }
  
  inputRDT <- file.path("data", paste(matrixName, "rarefaction-table.csv", sep="-"))
  if(verbose) 
    cat("\nUpload rarefaction table : from", inputRDT, "\n") 
  
  rarefactionTable <- read.csv(file=inputRDT, head=TRUE, sep=",", row.names=paste(levels, qs, sep=""), check.names=FALSE)
}

###### dissimilarity matrix #####
# Dissimilarity matrix of paired samples
# diss.fun = "beta1-1", "jaccard", "horn.morisita"
getDissimilarityMatrix <- function(expId, isPlot, min2, diss.fun="beta1-1", taxa.group="all", verbose=TRUE) {
  n <- length(matrixNames) 
  # hard code for Vegetation that only has plot and always keep singletons
  if (expId == n) {
    fname <- paste(matrixNames[expId], postfix("all", TRUE, FALSE, sep="-"), diss.fun, sep = "-")
  } else {
    fname <- paste(matrixNames[expId], postfix(taxa.group, isPlot, min2, sep="-"), diss.fun, sep = "-") 
  }
  
  inputB <- file.path("data", "dist", paste(fname, "csv", sep = "."))
  if(verbose) 
    cat("\nUpload", diss.fun, "matrix of", taxa.group, "taxa group(s) from", inputB, "\n") 
  
  diss.matrix <- readFile(file=inputB, sep=",")
  
  return(diss.matrix)
}

###### table to max remained diversity ##### 
getMaxRemainedDiversity <- function(lev.q, taxa.group="assigned", verbose=TRUE) {
  inputT <- file.path("data", "maxrd", paste("max-div", lev.q, taxa.group,"table.csv", sep = "-"))
  if (file.exists(inputT)) {
    max.rd <- read.csv(file=inputT, head=TRUE, sep=",", row.names=1, check.names=FALSE)
    if(verbose) 
      cat("\nUpload max remained diversity table from", inputT, "\n") 
  } else {
    max.rd <- NULL
    cat("Warning: cannot find max remained diversity table", inputT, "\n") 
  }
  max.rd
}

getMaxRemainedDiversityRank <- function(lev.q, taxa.group="assigned", verbose=TRUE) {
  inputT <- file.path("data", "maxrd", paste("max-div-rank", lev.q, taxa.group,"table.csv", sep = "-"))
  if (file.exists(inputT)) {
    max.rd <- read.csv(file=inputT, head=TRUE, sep=",", row.names=1, check.names=FALSE)
    if(verbose) 
      cat("\nUpload max remained diversity rank table from", inputT, "\n") 
  } else {
    max.rd <- NULL
    cat("Warning: cannot find max remained diversity rank table", inputT, "\n") 
  }
  max.rd
}

#' meta data of samples
#' 
#' @return \code{getSampleMetaData} returns a data frame
#' @export
#' @examples 
#' env <- getSampleMetaData(isPlot=TRUE)
#' 
#' @rdname getData
getSampleMetaData <- function(isPlot, verbose=TRUE) {
  if (isPlot) {
    inputCM <- file.path("data", "env", "LBI_all_env_data_by_plot.txt")
  } else {
    # e.g. data/16S.txt
    inputCM <- file.path("data", "env", "LBI_all_env_data_by_subplot.txt")
  }
  if(verbose) 
    cat("\nUpload enviornmental data from", inputCM, "\n") 
  env <- readFile(inputCM)
  
  env[,"ForestType"] <- gsub(":.*", "", env[,"ForestType"], ignore.case = T)
  env[,"ForestType"] <- gsub("x", "unknown", env[,"ForestType"], ignore.case = T)
}


######## elevations #######
getElevPlotDist <- function(plot.names, env.byplot) { 
  colElev = 1
  # case insensitive
  matched.id <- match(tolower(plot.names), tolower(rownames(env.byplot)))
  matched.id <- matched.id[!is.na(matched.id)]
  # match 
  env.plot.match <- env.plot[matched.id, ]
  
  cat("Find", nrow(env.plot.match), "plots having elevations, community matrix has", 
      length(plot.names), "plots, meta-data file has", nrow(env.plot), "plots.\n")
  
  return(dist(env.plot.match[,colElev]))
}


