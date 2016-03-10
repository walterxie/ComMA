# Picante: R tools for integrating phylogenies and ecology

#' @name CommunityPhyloStru
#' @title Community Phylogenetic Structure from \pkg{picante} Package
#'
#' @description Data input \strong{t.communityMatrix} is 
#' a transposed matrix from community matrix we defined in \pkg{ComMA}.
#' \strong{phyloTree} is a rooted tree of phylo object, 
#' which can get from \pkg{ape} \code{\link{read.tree}}.
#' 
#' \code{phylo.alpha} returns \pkg{picante} \code{\link{pd}}, 
#' which is phylogenetic alpha diversity (PD) index proposed by Faith (1992).
#' @param t.communityMatrix A transposed matrix from community matrix, where rows are samples, columns are OTUs
#' @param phyloTree A rooted tree of phylo object
#' @param ORD.RES The function how to order the sample names in the result
#' @param verbose default TRUE
#' @return 
#' \code{phylo.alpha} returns a data frame  of results from \code{\link{pd}}. 
#' @export
#' @examples 
#' pd.alpha <- phylo.alpha(t.communityMatrix, phyloTree)
#' 
#' @rdname CommunityPhyloStru
phylo.alpha <- function(t.communityMatrix, phyloTree, ORD.RES=function(res) {res[order(rownames(res)),]}, verbose=TRUE) {
  if ( ! all( sort(colnames(t.communityMatrix)) == sort(phyloTree$tip.label) ) ) 
    stop( paste("Community OTU names do not match tree tip labels") )
  
  if(verbose) {
    cat("Input community", nrow(t.communityMatrix), "samples", ncol(t.communityMatrix), "OTUs", 
        ", phylogenetic tree with", length(phyloTree$tip.label), "tips and", phyloTree$Nnode, "internal nodes.\n") 
    cat("Analysis: Faith's phylogenetic alpha diversity.\n") 
  }
  
  require(picante)
  
  # phylogenetic alpha diversity (PD) index proposed by Faith (1992)
  pd.result <- pd(t.communityMatrix, phyloTree, include.root = TRUE)
  
  pd.result <- ORD.RES(pd.result)
  pd.result
}

#' \code{phylo.mpd} returns \pkg{picante} \code{\link{ses.mpd}}, 
#' which is MPD standardized effect size of mean pairwise distances in communities.
#' When used with a phylogenetic distance matrix, equivalent to -1 times the Nearest Taxon Index (NTI).
#' @return 
#' \code{phylo.mpd} returns a data frame of results from \code{\link{ses.mpd}}. 
#' @export
#' @examples 
#' pd.mpd <- phylo.mpd(t.communityMatrix, phyloTree)
#' 
#' @rdname CommunityPhyloStru
phylo.mpd <- function(t.communityMatrix, phyloTree, ORD.RES=function(res) {res[order(rownames(res)),]}, verbose=TRUE) {
  if(verbose) 
    cat("Analysis: mean pairwise distance (MPD).\n")
  
  require(picante)
  
  # cophenetic distances for a hierarchical clustering
  phydist <- cophenetic(phyloTree)
  
  # When used with a phylogenetic distance matrix, equivalent to -1 times the Nearest Taxon Index (NTI).
  # MPD: standardized effect size of mean pairwise distances in communities
  ses.mpd.result <- ses.mpd(t.communityMatrix, phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 99)
  ses.mpd.result <- ORD.RES(ses.mpd.result)
  ses.mpd.result
}

#' \code{phylo.mntd} returns \pkg{picante} \code{\link{ses.mntd}}, 
#' which is MNTD standardized effect size of mean nearest taxon distances in communities
#' @return 
#' \code{phylo.mntd} returns a data frame of results from \code{\link{ses.mntd}}. 
#' @export
#' @examples 
#' pd.mntd <- phylo.mntd(t.communityMatrix, phyloTree)
#' @rdname CommunityPhyloStru
phylo.mntd <- function(t.communityMatrix, phyloTree, ORD.RES=function(res) {res[order(rownames(res)),]}, verbose=TRUE) {
  if(verbose) 
    cat("Analysis: mean nearest taxon distance (MNTD).\n")

  require(picante)
  
  # cophenetic distances for a hierarchical clustering
  phydist <- cophenetic(phyloTree)
  # MNTD: standardized effect size of mean nearest taxon distances in communities
  ses.mntd.result <- ses.mntd(t.communityMatrix, phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 99)
  ses.mntd.result <- ORD.RES(ses.mntd.result)
  ses.mntd.result
}

#' \code{phylo.beta.dist} returns dist object from \pkg{picante} \code{\link{comdist}}, 
#' which is phylogenetic beta diversity (Steven Kembel).
#' @return 
#' \code{phylo.beta.dist} returns a \code{\link{dist}} of results from \code{\link{comdist}}. 
#' @export
#' @examples 
#' pd.beta.dist <- phylo.beta.dist(t.communityMatrix, phyloTree)
#' 
#' @rdname CommunityPhyloStru
phylo.beta.dist <- function(t.communityMatrix, phyloTree, verbose=TRUE) {
  if(verbose) 
    cat("Analysis: MPD phylogenetic beta diversity.\n")
  
  require(picante)

  # cophenetic distances for a hierarchical clustering
  phydist <- cophenetic(phyloTree)
  # MPD (mean pairwise distance) separating taxa in two communities, phylogenetic beta diversity (Steven Kembel)
  comdist.result <- comdist(t.communityMatrix, phydist)
  comdist.result
}

#' \code{printResult} prints intermediate data to either file or console.
#' @param pd.alpha The result of phylogenetic alpha diversity, ignore it if NULL. 
#' @param pd.mpd The result of MPD, ignore it if NULL. 
#' @param pd.mntd The result of MNTD, ignore it if NULL. 
#' @param pd.beta.dist The result of phylogenetic beta diversity, ignore it if NULL.
#' @param filePath The path to create the file if tableFile is not NULL.
#' @param fileStem The file prefix to identify the result multi-conditions from multi-datasets
#' @param tableFile If NULL, then print the results to console, otherwise print them to the file. Default to NULL. 
#' @export
#' @examples 
#' tableFile <- file.path(workingPath, "report.tex")
#' pdFilePath <- file.path(workingPath, "data", "pd")
#' mkdir(pdFilePath) 
#' fileStem <- paste("16S", postfix(taxa.group, isPlot, sep="-"), sep = "-")
#' printResult(pd.alpha, pd.mpd, pd.mntd, pd.beta.dist, pdFilePath, fileStem, tableFile)
#' 
#' @rdname CommunityPhyloStru
printResult <- function(pd.alpha=NULL, pd.mpd=NULL, pd.mntd=NULL, pd.beta.dist=NULL, filePath, fileStem, tableFile=NULL) {
  if (!is.null(pd.alpha)) {
    fname <- file.path(filePath, paste(fileStem, "pd-alpha", sep = "-"))
    writeTable(pd.alpha, paste(fname, "csv", sep = "."))
    printXTable(pd.alpha, 
                caption = paste("Faith's phylogenetic alpha diversity of", fileStem), 
                label = paste("tab:pd:alpha", fileStem, sep = ":"), file=tableFile)    
  }
  
  if (!is.null(pd.mpd)) {
    fname <- file.path(filePath, paste(fileStem, "mpd", sep = "-"))
    writeTable(pd.mpd, paste(fname, "csv", sep = "."))
    printXTable(pd.mpd, 
                caption = paste("The mean pairwise distance (MPD) between all species in each community of", fileStem), 
                label = paste("tab:mpd:", fileStem, sep = ":"), file=tableFile)    
  }
  
  if (!is.null(pd.mntd)) {
    fname <- file.path(filePath, paste(fileStem, "mntd", sep = "-"))
    writeTable(pd.mntd, paste(fname, "csv", sep = "."))
    printXTable(pd.mntd, 
                caption = paste("The mean nearest taxon distance (MNTD)", 
                                "separating each species in the community from its closest relative of", fileStem), 
                label = paste("tab:mntd:", fileStem, sep = ":"), file=tableFile) 
  }
  
  if (!is.null(pd.beta.dist)) {
    comdist.m <- as.matrix(pd.beta.dist)
    
    fname <- file.path(filePath, paste(fileStem, "pd-beta", sep = "-"))
    writeTable(comdist.m, paste(fname, "csv", sep = "."))
    
    comdist.m <- comdist.m[order(rownames(comdist.m)),]
    comdist.m <- comdist.m[,order(colnames(comdist.m))]
    comdist.m[upper.tri(comdist.m)] <- NA
    comdist.m[comdist.m=="0"] <- NA
    comdist.m <- comdist.m[-1,-ncol(comdist.m)]
    
    printXTable(comdist.m, 
                caption = paste("Phylogenetic beta diversity of", fileStem), 
                label = paste("tab:pd:beta", fileStem, sep = ":"), file=tableFile) 
  }
}

#' \code{plotPDBeta} plots the dendrogram of phylogenetic beta diversity (Steven Kembel)
#' @export
#' @examples 
#' figDir <- file.path(workingPath, "figures")
#' fileStem <- paste("16S", postfix(taxa.group, isPlot, sep="-"), sep = "-")
#' plotPDBeta(pd.beta.dist, figDir, fileStem)
#' 
#' @rdname CommunityPhyloStru
plotPDBeta <- function(pd.beta.dist, filePath, fileStem) {
  comdist.clusters <- hclust(pd.beta.dist)
  
  fname <- file.path(filePath, fileStem)
  pdf(paste(fname, "pdf", sep = "."), width=10, height=5)
  plot(comdist.clusters, xlab="", sub ="", main=paste("Phylogenetic beta diversity"))
  invisible(dev.off())
}


