
#' @name getData
#' @title Get example data for \pkg{ComMA} package
#'
#' @description The functions are used to setup data for the project. 
#' The data format is described in \code{\link{ComMA}}.
#' 
#' @details 
#' \code{postfix} adds postfix for figure name or table label, 
#' and returns a name concatenated by substrings and separated by \code{sep}.
#' @param ... A list of names to be concatenated.
#' @param sep Default to dash "-".
#' @export
#' @examples 
#' postfix("16S", "assigned")
#'[1] "16S-assigned-byplot-min2"
#'
#' @rdname getData
postfix <- function(..., isPlot=TRUE, taxa.group="all", minAbund=2, minRich=1, sep="-") {
  full.name <- paste(list(...), collapse = sep)
  
  if (taxa.group!="all") 
    full.name <- paste(full.name, taxa.group, sep = sep) 
  if (isPlot) 
    full.name <- paste(full.name, "byplot", sep = sep) 
  if (minAbund > 1) 
    full.name <- paste(full.name, paste0("min", minAbund), sep = sep)
  if (minRich > 1) 
    full.name <- paste(full.name, minRich, sep = sep)
  
  return(full.name)
}

#' @details \code{getPlot} extracts plot names from full names (plot + subplot) separated by \code{sep}.
#' 
#' @param full.name The full name has plot and subplot together, but separated by \code{sep}.
#' @export
#' @examples 
#' getPlot(c("Plot1-A", "CM30c39-L"))
#' [1] "Plot1" "CM30c39"
#' 
#' @rdname getData
getPlot <- function(full.name, sep="-") {
  sapply(strsplit(as.character(full.name), sep), "[[", 1)
}

# Dropbox: Hauturu Miseq analysis data
######## load community matrix #######
#' @details \code{getCommunityMatrix} returns a community matrix, 
#' where rows are OTUs or individual species and columns are sites or samples. 
#' 
#' @param matrix.name The string to locate the matrix from its file name.
#' @param isPlot Boolean value to determine the matrix file sampled by subplot or plot
#' @param minAbund The minimum abundance threshold to remove rows/columns 
#' by row/column sum of abundance. For exampe, if minAbund=2, then remove 
#' all singletons appeared in only one sample. If minAbund=1, 
#' then remove all empty rows/columns. Default to 2 (singletons).
#' But \code{postfix} is only used for naming, no data processed.
#' @param verbose More details. Default to TRUE.
#' @export
#' @examples 
#' # keep singletons
#' community.matrix <- getCommunityMatrix("16S", isPlot=TRUE, minAbund=1)
#' 
#' @rdname getData
getCommunityMatrix <- function(matrix.name, isPlot, minAbund=2, verbose=TRUE) {
  if (isPlot) {
    inputCM <- file.path("data", "OTU_tables", paste(matrix.name, "otutable_by_plot.txt", sep="_"))
  } else {
    # e.g. data/16S.txt
    inputCM <- file.path("data", "OTU_tables", paste(matrix.name, "otutable.txt", sep=""))
  }
  
  community.matrix <- ComMA::readCommunityMatrix(inputCM, matrix.name = matrix.name, 
                                                minAbund = minAbund, verbose = verbose)
  
  return(community.matrix)
}

#' @details 
#' \code{getCommunityMatrixT} returns a transposed matrix of (maybe also a subset of) 
#' community matrix, given more filters, such as \code{taxa.group}, \code{minRich}.
#' Its columns are OTUs or individual species and rows are sites or samples. 
#' It is also the abundances argument in \pkg{vegetarian} \code{\link{d}}.
#' return \emph{NULL}, if no OTUs or OTUs less than \code{minRich} threshold.
#' 
#' @param minRich The minimum richness to keep matrix. For example, 
#' drop the matrix (return NULL) if OTUs less than this threshold. Default to 1.
#' @export
#' @examples 
#' # by plot, remove singletons, BACTERIA only
#' t.community.matrix <- getCommunityMatrixT("16S", isPlot=TRUE, minAbund=2, taxa.group="BACTERIA")
#' 
#' @rdname getData
getCommunityMatrixT <- function(matrix.name, isPlot, taxa.group="all", minAbund=2, minRich=1, verbose=TRUE) {
  community.matrix <- ComMA::getCommunityMatrix(matrix.name, isPlot, minAbund)
  
  if (taxa.group != "all") {
    ##### load data #####
    taxaPaths <- ComMA::getTaxaPaths(matrix.name, taxa.group)
    
    if (nrow(taxaPaths) < minRich) {
      warning("Return NULL, because taxa table has ", nrow(taxaPaths), " row(s) classified as ", 
              taxa.group, " < ", minRich, " threshold.\n")
      return(NULL)
    } else {
      # merge needs at least 2 cols 
      taxaAssgReads <- merge(community.matrix, taxaPaths, by = "row.names")
      
      if (nrow(taxaAssgReads) < minRich) {
        warning("Return NULL, because cm has ", nrow(taxaAssgReads), " row(s) match ", taxa.group, 
                " < ", minRich, " threshold.\n")
        return(NULL)
      }
      
      # move 1st col Row.names to row.names
      rownames(taxaAssgReads) <- taxaAssgReads[,"Row.names"]
      taxaAssgReads <- taxaAssgReads[,-1]
      # get CM
      taxaAssgReads <- taxaAssgReads[,1:ncol(community.matrix)]
      
      cat("Merging", nrow(taxaAssgReads), "matched OTUs from", nrow(community.matrix), "OTUs in matrix to", 
          nrow(taxaPaths), "taxa classification, taxa.group =", taxa.group, ", final ncol =", ncol(taxaAssgReads), ".\n")
      
      community.matrix <- data.matrix(taxaAssgReads)
    }
  }
  t.community.matrix <- ComMA::transposeCM(community.matrix)
  
  return(t.community.matrix)
}
# rownames(community.matrix) <- gsub("-(.*)|", "\\1", rownames(community.matrix))

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
#' @details 
#' \code{getTaxaPaths} returns a taxonomic matrix, 
#' where rows are OTUs or individual species (they have to be the subset of rows 
#' in the community matrix from \code{getCommunityMatrix}), 
#' and columns are taxonomic ranks, such as c("superkingdom", "kingdom", "phylum", 
#' "class", "order", "family", "genus", "species") inlcuding full taxonomic path. 
#' @export
#' @examples 
#' taxaPaths <- getTaxaPaths(matrix.name="16S", taxa.group="BACTERIA|ARCHAEA")
#' taxaPaths <- getTaxaPaths(matrix.name="16S", taxa.group="BACTERIA|ARCHAEA", include=FALSE)
#' taxaPaths <- getTaxaPaths(matrix.name="18S", taxa.group="PROKARYOTA", rank="superkingdom")
#' 
#' @rdname getData
getTaxaPaths <- function(matrix.name, taxa.group="all", rank="kingdom", 
                         include=TRUE, regex="(\\|[0-9]+)", verbose=TRUE) {
  inputTaxa <- file.path("data", "Taxonomy_tables", paste(matrix.name, "taxonomy_table.txt", sep="_"))
  taxaPaths <- ComMA::readTaxaTable(inputTaxa, matrix.name=matrix.name, taxa.group=taxa.group, 
                                    rank=rank, include=include, regex=regex, verbose=verbose)
  
  nTaxa=nrow(taxaPaths)
  ##### keep OTU rows contain given taxa belongTo ##### 
  if (taxa.group != "all") {
    # Exclude probably bogus taxa
    taxaPaths <- subset(taxaPaths, !(grepl("Cnidaria|Brachiopoda|Echinodermata|Porifera", taxaPaths[,"phylum"])|
                                       grepl("Bivalvia|Teleostei|Elasmobranchii|Polyplacophora", taxaPaths[,"class"])|
grepl("Nudibranchia|Crocodylia|Serpentes|Testudines|Carnivora|Gymnophiona|Lagomorpha|Rodentia|Serpentes|Scorpiones", taxaPaths[,"order"])))
  }
  
  if (nrow(taxaPaths) < 1)
    warning("Cannot find ", taxa.group, " from taxa path file ", inputTaxa, " !")
  
  if(verbose && nrow(taxaPaths) < nTaxa) 
    cat("\nSelect", nrow(taxaPaths), "taxa classification after excluding probably bogus taxa.\n") 
  
  return(taxaPaths)
}


# rankLevel: the taxa level in each bar
# groupLevel: used to assign colour for each group, and must higher than rankLevel
# taxa.group: keep OTU rows contain given taxa group, if "all", keep all
# return taxaAssgReads = CM + rankLevel + groupLevel
getTaxaAssgReads <- function(matrix.name, isPlot, minAbund=2, rankLevel, groupLevel, taxa.group="all") {
  cat("Create taxonomy assignment for", matrix.name, ".\n")
  
  ##### load data #####
  community.matrix <- ComMA::getCommunityMatrix(matrix.name, isPlot, minAbund)
  
  community.matrix <- community.matrix[order(rownames(community.matrix)),]
  community.matrix <- community.matrix[,order(colnames(community.matrix))]
  
  taxaPaths <- ComMA::getTaxaPaths(matrix.name, taxa.group)
  
  ###### taxa assignment by reads #####
  if ( ! tolower(rankLevel) %in% tolower(colnames(taxaPaths)) ) 
    stop( paste("Column name", rankLevel, "not exist in taxa path file for ", matrix.name) )
  if (! tolower(groupLevel) %in% tolower(colnames(taxaPaths)) ) 
    stop( paste("Column name", groupLevel, "not exist in taxa path file for ", matrix.name) )
  
  colRankLevel <- which(tolower(colnames(taxaPaths))==tolower(rankLevel))
  colGroupLevel <- which(tolower(colnames(taxaPaths))==tolower(groupLevel))
  
  taxaAssgReads <- merge(community.matrix, taxaPaths[,c(colRankLevel, colGroupLevel)], by = "row.names")
  
  taxaAssgReads[,rankLevel] <- gsub("(\\s\\[=.*\\])", "", taxaAssgReads[,rankLevel])
  taxaAssgReads[,groupLevel] <- gsub("(\\s\\[=.*\\])", "", taxaAssgReads[,groupLevel])
  
  cat("Merging:", nrow(taxaAssgReads), "OTUs are matched from", nrow(community.matrix), "OTUs in matrix to", 
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

#' @details \code{getSampleMetaData} returns a data frame containing meta data of samples
#' 
#' @export
#' @examples 
#' env <- getSampleMetaData(isPlot=TRUE)
#' 
#' @rdname getData
getSampleMetaData <- function(isPlot, verbose=TRUE) {
  if (isPlot) {
    inputCM <- file.path("data", "Environmental_data", "LBI_all_env_data_by_plot.txt")
  } else {
    # e.g. data/16S.txt
    inputCM <- file.path("data", "Environmental_data", "LBI_all_env_data_by_subplot.txt")
  }
  env <- ComMA::readFile(inputCM, verbose=verbose, msg.file="enviornmental data", msg.row="samples")
  
  env[,"ForestType"] <- gsub(":.*", "", env[,"ForestType"], ignore.case = T)
  env[,"ForestType"] <- gsub("x", "unknown", env[,"ForestType"], ignore.case = T)
  
  return(env)
}

#' @details \code{getPhyloTree} returns a rooted tree of \code{\link{phylo}} object.
#' 
#' @export
#' @examples 
#' phylo.tree <- getPhyloTree("16S", "bacteria")
#' 
#' @rdname getData
getPhyloTree <- function(matrix.name, taxa.group="assigned", minAbund=2, verbose=TRUE) {
  inputT <- file.path("data","Trees",paste(postfix(matrix.name, taxa.group=tolower(taxa.group), 
                                                   minAbund=minAbund, isPlot=FALSE), "tre", sep = "."))
  if (file.exists(inputT)) {
    require(ape)
    cat("Load tree from", inputT, "\n") 
    tree <- read.tree(inputT)
    if(verbose) print(tree)
  } else {
    tree <- NULL
    warning("Cannot find tree file: ", inputT, "\n") 
  }
  tree
}

###### Intermediate Data ##### 

#' table to plot Phylo Rarefaction
getPhyloRareTable <- function(expId, isPlot, min2, taxa.group="assigned", verbose=TRUE) {
  n <- length(matrixNames) 
  # hard code for Vegetation that only has plot and always keep singletons
  if (expId == n) {
    mid.name <- ComMA::postfix("all", TRUE, FALSE, sep="-")
  } else {
    mid.name <- ComMA::postfix(taxa.group, isPlot, min2, sep="-") 
  }
  
  inputT <- file.path("data", "pdrf", paste(matrixNames[expId], mid.name, "phylorare", "table.csv", sep="-"))
  if (file.exists(inputT)) {
    phylo.rare.df <- read.csv(file=inputT, head=TRUE, sep=",", row.names=1, check.names=FALSE)
    if(verbose) 
      cat("\nUpload phylo rarefaction table from", inputT, "\n") 
  } else {
    phylo.rare.df <- NULL
    warning("Cannot find phylo rarefaction table ", inputT, " \n") 
  }
  phylo.rare.df
}

#' table to plot Rarefaction
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
    warning("Cannot find rarefaction table per sample ", inputT, " \n") 
  }
  rare.df
}

getRarefactionTable <- function(expId, isPlot, min2, verbose=TRUE) {
  n <- length(matrixNames) 
  matrixName <- matrixNames[expId]
  # hard code for Vegetation that only has plot and always keep singletons
  if (expId == n) {
    matrixName <- ComMA::postfix(matrixName, TRUE, FALSE, sep="-")
  } else {
    matrixName <- ComMA::postfix(matrixName, isPlot, min2, sep="-") 
  }
  
  inputRDT <- file.path("data", paste(matrixName, "rarefaction-table.csv", sep="-"))
  if(verbose) 
    cat("\nUpload rarefaction table : from", inputRDT, "\n") 
  
  rarefactionTable <- read.csv(file=inputRDT, head=TRUE, sep=",", row.names=paste(levels, qs, sep=""), check.names=FALSE)
}

#' dissimilarity matrix
#' Dissimilarity matrix of paired samples
#' diss.fun = "beta1-1", "jaccard", "horn.morisita"
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
  
  diss.matrix <- ComMA::readFile(file=inputB, sep=",")
  
  return(diss.matrix)
}

#' table to max remained diversity
getMaxRemainedDiversity <- function(lev.q, taxa.group="assigned", verbose=TRUE) {
  inputT <- file.path("data", "maxrd", paste("max-div", lev.q, taxa.group,"table.csv", sep = "-"))
  if (file.exists(inputT)) {
    max.rd <- read.csv(file=inputT, head=TRUE, sep=",", row.names=1, check.names=FALSE)
    if(verbose) 
      cat("\nUpload max remained diversity table from", inputT, "\n") 
  } else {
    max.rd <- NULL
    warning("Cannot find max remained diversity table ", inputT, " \n") 
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
    warning("Cannot find max remained diversity rank table ", inputT, " \n") 
  }
  max.rd
}

#' elevations
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


