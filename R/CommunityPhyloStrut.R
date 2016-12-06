# Picante: R tools for integrating phylogenies and ecology

#' @name CommunityPhyloStru
#' @title Community Phylogenetic Structure from \pkg{picante} Package
#'
#' @description A wrapper class using \pkg{picante} Package. 
#' Data input \strong{t.community.matrix} is 
#' a transposed matrix from community matrix we defined in \pkg{ComMA}.
#' \strong{phylo.tree} is a rooted tree of phylo object, 
#' which can get from \pkg{ape} \code{\link{read.tree}}.
#' 
#' @details \code{phylo.alpha} returns a data frame of the PD and species richness (SR) 
#' values for all samples from \code{\link{pd}} in \pkg{picante}. 
#' Phylogenetic alpha diversity (PD) index is proposed by Faith (1992).
#' @param t.community.matrix A transposed matrix from community matrix, 
#' where rows are samples, columns are OTUs.
#' @param phylo.tree A rooted tree of phylo object
#' @param ORD.RES The function how to order the sample names in the result
#' @param verbose default TRUE
#' @note If taxa in phylogenies do not match OTUs in the community, 
#' then use \code{\link{match.phylo.comm}} before \code{phylo.alpha}. 
#' See examples. But prefer \code{comm > phy} in practical, 
#' otherwise \code{combined$phy} will be unrooted tree. 
#' Therefore, in the case of single-species samples the PD will 
#' be equal to NA (include.root=FALSE), see \code{\link{pd}}.
#' @export
#' @examples 
#' pd.alpha <- phylo.alpha(t.community.matrix, phylo.tree)
#' pd.alpha
#' \tabular{rrr}{
#'    \tab PD \tab SR\cr
#'   CM30b51 \tab 402.0822 \tab 2796\cr
#'   CM30b58 \tab 456.4397 \tab 2634\cr
#'   CM30b60 \tab 554.9873 \tab 2931
#' }
#' 
#' combined <- match.phylo.comm(phylo.tree, t.communityMatrix)
#' pd.alpha <- phylo.alpha(combined$comm, combined$phy, ...) 
#' 
#' @rdname CommunityPhyloStru
phylo.alpha <- function(t.community.matrix, phylo.tree, include.root = TRUE,
                        ORD.RES=function(res) {res[order(rownames(res)),]}, verbose=TRUE) {
  if(verbose) {
    cat("Input community", nrow(t.community.matrix), "samples", ncol(t.community.matrix), "OTUs", 
        ", phylogenetic tree with", length(phylo.tree$tip.label), "tips and", phylo.tree$Nnode, "internal nodes.\n") 
    cat("Analysis: Faith's phylogenetic alpha diversity.\n") 
  }
  if ( ! all( sort(colnames(t.community.matrix)) == sort(phylo.tree$tip.label) ) ) 
    warning("community OTU names do not match tree tip labels !",
            "\nPlease apply match.phylo.comm before this.")
  
  require(picante)
  if (!is.rooted(phylo.tree) && include.root) {
    warning("phylo.tree is unrooted, set include.root=FALSE !")
    include.root=FALSE
  }
  # phylogenetic alpha diversity (PD) index proposed by Faith (1992)
  pd.result <- pd(t.community.matrix, phylo.tree, include.root = include.root)
  pd.result <- ORD.RES(pd.result)
  pd.result
}

#' @details \code{phylo.mpd} returns \pkg{picante} \code{\link{ses.mpd}}, 
#' which is MPD standardized effect size of mean pairwise distances in communities.
#' When used with a phylogenetic distance matrix, equivalent to -1 times the Nearest Taxon Index (NTI).
#' @export
#' @examples 
#' pd.mpd <- phylo.mpd(t.community.matrix, phylo.tree)
#' 
#' @rdname CommunityPhyloStru
phylo.mpd <- function(t.community.matrix, phylo.tree, ORD.RES=function(res) {res[order(rownames(res)),]}, verbose=TRUE) {
  if(verbose) 
    cat("Analysis: mean pairwise distance (MPD).\n")
  
  require(picante)
  
  # cophenetic distances for a hierarchical clustering
  phydist <- cophenetic(phylo.tree)
  
  # When used with a phylogenetic distance matrix, equivalent to -1 times the Nearest Taxon Index (NTI).
  # MPD: standardized effect size of mean pairwise distances in communities
  ses.mpd.result <- ses.mpd(t.community.matrix, phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 99)
  ses.mpd.result <- ORD.RES(ses.mpd.result)
  ses.mpd.result
}

#' @details \code{phylo.mntd} returns \pkg{picante} \code{\link{ses.mntd}}, 
#' which is MNTD standardized effect size of mean nearest taxon distances in communities.
#' @export
#' @examples 
#' pd.mntd <- phylo.mntd(t.community.matrix, phylo.tree)
#' @rdname CommunityPhyloStru
phylo.mntd <- function(t.community.matrix, phylo.tree, ORD.RES=function(res) {res[order(rownames(res)),]}, verbose=TRUE) {
  if(verbose) 
    cat("Analysis: mean nearest taxon distance (MNTD).\n")

  require(picante)
  
  # cophenetic distances for a hierarchical clustering
  phydist <- cophenetic(phylo.tree)
  # MNTD: standardized effect size of mean nearest taxon distances in communities
  ses.mntd.result <- ses.mntd(t.community.matrix, phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 99)
  ses.mntd.result <- ORD.RES(ses.mntd.result)
  ses.mntd.result
}

#' @details \code{phylo.beta.dist} returns \code{\link{dist}} object from \pkg{picante} \code{\link{comdist}}, 
#' which is phylogenetic beta diversity (Steven Kembel).
#' @export
#' @examples 
#' pd.beta.dist <- phylo.beta.dist(t.community.matrix, phylo.tree)
#' 
#' @rdname CommunityPhyloStru
phylo.beta.dist <- function(t.community.matrix, phylo.tree, verbose=TRUE) {
  if(verbose) {
    cat("Analysis: MPD phylogenetic beta diversity.\n")
    cat("t.community.matrix has", nrow(t.community.matrix), "samples,", 
        ncol(t.community.matrix), "OTUs/taxa,", ".\n")
    cat("phylo.tree has ", length(phylo.tree$tip.label), " tips and ", 
        phylo.tree$Nnode, " internal nodes.\n")
  }
  require(picante)

  # cophenetic distances for a hierarchical clustering
  phydist <- cophenetic(phylo.tree)
  # MPD (mean pairwise distance) separating taxa in two communities, phylogenetic beta diversity (Steven Kembel)
  comdist.result <- comdist(t.community.matrix, phydist)
  comdist.result
}

#' @details \code{phylo.unifrac.dist} calculates unweighted UniFrac and 
#' returns dist object from \pkg{picante} \code{\link{unifrac}}.
#' @export
#' @examples 
#' unifrac.dist <- phylo.unifrac.dist(t.community.matrix, phylo.tree)
#' 
#' @rdname CommunityPhyloStru
phylo.unifrac.dist <- function(t.community.matrix, phylo.tree, verbose=TRUE) {
  if(verbose) {
    cat("Analysis: calculates unweighted UniFrac.\n")
    cat("t.community.matrix has", nrow(t.community.matrix), "samples,", 
        ncol(t.community.matrix), "OTUs/taxa,", ".\n")
    cat("phylo.tree has ", length(phylo.tree$tip.label), " tips and ", 
        phylo.tree$Nnode, " internal nodes.\n")
  }
  require(picante)
  
  unifrac.dist <- unifrac(t.community.matrix, phylo.tree)
  unifrac.dist
}


#' @details \code{printResult} prints intermediate data to either file or console.
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

#' @details \code{plotPDBeta} plots the dendrogram of phylogenetic beta diversity (Steven Kembel)
#' @export
#' @examples 
#' plotPDBeta(pd.beta.dist)
#' pdf.plot(plot, fig.path="pd-beta.pdf", width=10, height=5)  
#' 
#' @rdname CommunityPhyloStru
plotPDBeta <- function(pd.beta.dist, title="Phylogenetic beta diversity", xlab="", sub ="", ...) {
  comdist.clusters <- hclust(pd.beta.dist)
  p <- plot(comdist.clusters, xlab=xlab, sub=sub, main=title, ...)
  return(p)
}
