# Author: Andrew Dopheide, Walter Xie
# Accessed on 25 Nov 2016


#' @name CommunityComparison
#' @title Compare community structure
#'
#' @description
#' Multivariate patterns of sample similarity were compared between pairwised community matrices 
#' using Procrustes and Mantel tests.
#' 
#' @details
#' \code{getDissimilarityList} calculates a list of similarity/dissimilarity 
#' given a list of community matrices and the distance metric.  
#' 
#' @param cm.list A list of community matrices.
#' @param metric,... The distance metric, default to "jaccard",
#' and other parametes passed to \code{\link{getDissimilarity}}.
#' @keywords comparison
#' @export
#' @examples 
#' jaccard.dist <- getDissimilarityList(cm.list, metric="jaccard")
#' jaccard.dist$dist.list
#' 
#' @rdname CommunityComparison
getDissimilarityList <- function(cm.list, metric="jaccard", ...){
  dist.list <- list()
  
  for (i in 1:length(cm.list)) {
    cm.name <- names(cm.list)[i]
    if (is.null(cm.name))
      cm.name <- paste0("community",i)
    cat("Calculate", metric, "dissimilarity for", cm.name, ".\n")
    
    dist.cm <- ComMA::getDissimilarity(cm.list[[i]], method=metric, ...)
    dist.list[[cm.name]] <- dist.cm
  }
  
  list(dist.list=dist.list, metric=metric)
}

#' @details
#' \code{mantelComparison} makes Mantel tests for comparisons 
#' between the overall community structure.  
#' 
#' @param dist.list A list of dissimilarity in \code{\link{dist}} 
#' between pairwised community matrices, 
#' which can be calculated by \code{\link{getDissimilarityList}}. 
#' @param method,permutations Parameters used in \code{\link{mantel}} in \pkg{vegan}.
#' @export
#' @examples 
#' mantel <- mantelComparison(jaccard.dist$dist.list)
#' mantel.tri <- getTriMatrix(mantel$m.df) 
#' 
#' @rdname CommunityComparison
mantelComparison <- function(dist.list, method="pearson", permutations = 999){
  m.list <- list() # list for mantel results storage

  label.matched <- getLabelMatchedDist(dist.list)
  require(vegan)
  for(i in 1:length(label.matched$dist1)){
    man <- mantel(label.matched$dist1[[i]], label.matched$dist2[[i]], 
                  method=method, permutations = permutations)
    l1 = label.matched$pairs[[i]][[1]] 
    l2 = label.matched$pairs[[i]][[2]]
    m.list[[i]] <- list(l1 = l1, l2 = l2, corr = man$statistic, sign = man$signif)
    cat("Pair (", l1, ",", l2, ") correlation =", man$statistic, ", significance =", man$signif, "\n")
  }
  m.result <- do.call("rbind", lapply(m.list, data.frame))
  
  list(m.df=m.result, m.list=m.list, pairs=label.matched$pairs, samples=label.matched$samples, 
       same.samples=label.matched$same.samples)
}

#' @details
#' \code{procrustesComparison} makes Procrustes comparisons 
#' between the overall community structure.  
#' 
#' @param scale,symmetric,permutations Parameters used in 
#' \code{\link{procrustes}} and \code{\link{protest}}.
#' @export
#' @examples 
#' procrustes <- procrustesComparison(jaccard.dist$dist.list) 
#' prot.tri <- getTriMatrix(procrustes$prot.df, order.by=order.by) 
#' 
#' # Mantel test (lower triangle) and Procrustes test (upper triangle)
#' corrs <- combineTriMatrix(mantel.tri, prot.tri)
#' 
#' @rdname CommunityComparison
procrustesComparison <- function(dist.list, scale = TRUE, symmetric = TRUE, permutations = 999){  
  proc.list <- list() # list to store procrustes objects (for plotting)
  prot.list <- list()
  
  label.matched <- getLabelMatchedDist(dist.list)
  require(vegan)
  for(i in 1:length(label.matched$dist1)){
    # Generate MDS plots
    d1.mds <- metaMDS(label.matched$dist1[[i]])
    d2.mds <- metaMDS(label.matched$dist2[[i]])
    proc <- procrustes(d1.mds, d2.mds, scale = scale, symmetric = symmetric)
    prot <- protest(d1.mds, d2.mds, scale = scale, symmetric = symmetric, 
                    permutations = permutations)
    proc.list[[i]] <- proc
    l1 = label.matched$pairs[[i]][[1]]
    l2 = label.matched$pairs[[i]][[2]]
    prot.list[[i]] <- list(l1 = l1, l2 = l2, ss=prot$ss, corr=prot$t0, sign=prot$signif)
    cat("Pair (", l1, ",", l2, ") correlation =", prot$t0, ", significance =", prot$signif, 
        ", sum of squares ss =", prot$ss,"\n")
  }
  p.result <- do.call("rbind", lapply(prot.list, data.frame))
  
  list(proc.list = proc.list, prot.df=p.result, prot.list=prot.list, pairs=label.matched$pairs, 
       samples=label.matched$samples, same.samples=label.matched$same.samples)
}

#' @details
#' \code{sublistByPairs} takes a subset result from comparisons 
#' given the subset of pairwised communities.  
#' 
#' @param m.p.list Either \code{m.list} of the output from \code{mantelComparison}, 
#' or \code{proc.list} or \code{prot.list} from \code{procrustesComparison}.
#' @param pairs,subset.pairs The list of the list of two pairs of community names, 
#' see example. 
#' The string names have to be exaclty same in order to find matches, 
#' but the order of the pair does not matter the 1st or 2nd come first.
#' @export
#' @examples 
#' sub.proc <- sublistByPairs(procrustes$proc.list, procrustes$pairs, subset.pairs=list(list("16S bacteria","18S protists"),list("18S animals","COI−300 animals")) )
#' sub.proc$sub.list
#' 
#' @rdname CommunityComparison
sublistByPairs <- function(m.p.list, pairs, subset.pairs) {
  if (length(m.p.list) != length(pairs))
    stop("length(m.p.list) has to be same as length(pairs) !\n", 
         "m.p.list is either m.list of the output from function mantelComparison, ", 
         "or proc.list/prot.list from procrustesComparison.")
  m1 <- do.call("rbind", pairs)
  v1 <- paste(m1[,1], m1[,2])
  m2 <- do.call("rbind", subset.pairs)
  v2 <- paste(m2[,1], m2[,2])
  # check the reverse pair
  v3 <- paste(m2[,2], m2[,1])
  
  id.match2 <- match(tolower(v2), tolower(v1))
  id.match3 <- match(tolower(v3), tolower(v1))
  # merge two id vectors by non NA values
  id.match <- pmin(id.match2, id.match3, na.rm = TRUE)
  id.match <- id.match[!is.na(id.match)] 
  if (length(id.match) < 1) {
    warning("Cannot match subset.pairs to pairs !\n", 
            "subset.pairs = ", paste(subset.pairs, collapse = ","), 
            "\npais = ", paste(pairs, collapse = ","))
    return(list())
  }
  if (length(id.match) != length(v2))
    warning(length(id.match), " selected pairs != ", length(v2), " given subset pairs !\nSelected pairs are ", 
        paste(v1[id.match], collapse = ", "), ".\n")
  subset.list <- m.p.list[id.match]
  list(sub.list = subset.list, sub.pairs=subset.pairs)
}

#' @details
#' \code{plotProcrustes} plots Procrustes correlations between pairwised communities.  
#' 
#' @param proc.list The list of \code{\link{procrustes}} results.
#' @export
#' @examples 
#' plotProcrustes(procrustes$proc.list)
#' 
#' @rdname CommunityComparison
plotProcrustes <- function(proc.list, attr.df, colour.id="Elevation", 
                           title.list=list(), proc.list.pairs=list()) {
  if (length(title.list) > 0 && length(title.list) != length(proc.list))
    stop("length(title.list) has to be same as length(proc.list) !")
  if (length(proc.list.pairs) > 0 && length(proc.list.pairs) != length(proc.list))
    stop("length(proc.list.pairs) has to be same as length(proc.list) !")
  if (! colour.id %in% colnames(attr.df))
    stop("Invalid colour.id,", colour.id,  "not exsit in meta data column names !\n")
  
  require(ggplot2)
  plot.list <- list()
  for(i in 1:length(proc.list)){
    pro <- proc.list[[i]]
    pts <- data.frame(yMDS1 = pro$Yrot[,1], yMDS2 = pro$Yrot[,2], # rotated matrix i.e. d2.mds
                      xMDS1 = pro$X[,1], xMDS2 = pro$X[,2]) # target matrix i.e. d1.mds
    
    pts <- merge(pts, attr.df, by = "row.names")
    
    r1 = acos(pro$rotation[1,1]) # X axis rotation (radians)
    r2 = r1 + (pi/2) # Y axis rotation (radians)
    p <- ggplot(pts) +
      geom_point(aes_string(x = "yMDS1", y = "yMDS2", colour = colour.id), shape = 1.5, size = 1.5, alpha = 0.75) + # rotated i.e. d2 (circles)
      geom_point(aes_string(x = "xMDS1", y = "xMDS2", colour = colour.id), shape = 2, size = 1, alpha = 0.75) + # target i.e. d1 (triangles)
      scale_shape(solid = FALSE) + xlab("") + ylab("") +
      geom_segment(aes_string(x = "yMDS1", y = "yMDS2", xend = "xMDS1", yend = "xMDS2", colour = colour.id), alpha = 0.75) +
      #arrow = arrow(length = unit(0.2,"cm")), alpha = 0.75) + 
      #geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") + 
      #geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
      #geom_abline(intercept = 0, slope = tan(r1), colour = "grey") + 
      #geom_abline(intercept = 0, slope = tan(r2), colour = "grey") +
      #geom_text(aes(x = xMDS1, y = xMDS2, label = pts.mds$Row.names, colour = Elevation), size = 1.5, vjust = 1.5, alpha = 0.5) + 
      scale_colour_gradientn(colours = c("blue", "orange")) + labs(colour="Elevation (m)") +
      theme(panel.grid = element_blank(), plot.title = element_text(size = 8), 
            plot.margin = unit(c(0.1,0.1,0.1,0), "cm"), legend.key.width = unit(0.65, "cm")) 
   
    title <- NULL
    if (length(title.list) > 0) {
      title <- title.list[[i]]
    } else if (length(proc.list.pairs) > 0) {
      title <- paste0(letters[[i]], ". ", proc.list.pairs[[i]][[1]], " vs. ", proc.list.pairs[[i]][[2]]) #, ", ", metric, " distance")) #+
    } else if (!is.null(names(proc.list))) {     
      title <- paste0(letters[[i]], ". ", names(proc.list)[[i]])
    }
    if (!is.null(title))
      p <- p + ggtitle(title)
    
    #coord_fixed()
    tmp <- ggplot_gtable(ggplot_build(p))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    
    p <- p + theme(legend.position = "none")
    plot.list[[i]] <- p
    cat("Add Procrustes i=", i, ", title=", title, "\n")
  }
  return(list(plot.list, legend))
}



###### Internal #####

getLabelMatchedDist <- function(dist.list) {
  dist1.list <- list()
  dist2.list <- list()
  pair.list <- list()
  samples.list <- list()
  same.samples.list <- list()
  x <- combn(dist.list, 2, simplify = F) # Get all pairs of dist matrices
  for(i in 1:length(x)){
    cat("Paired dist matrices : ", labels(x[[i]]), "\n")
    if( !all(labels(x[[i]][[1]]) == labels(x[[i]][[2]])) ) { # Different numbers of samples
      cat("Making samples match...\n")
      # Subset to shared samples
      keep <- intersect(labels(x[[i]][[1]]), labels(x[[i]][[2]]))
      m1 <- as.matrix(x[[i]][[1]])
      m2 <- as.matrix(x[[i]][[2]])
      d1 <- as.dist(m1[keep, keep])
      d2 <- as.dist(m2[keep, keep])
      same.samples.list[[i]] <- FALSE
    } else { # Same numbers of samples
      d1 <- x[[i]][[1]]
      d2 <- x[[i]][[2]]
      same.samples.list[[i]] <- TRUE
    }
    cat("Samples match ... ", all(names(d1) == names(d2)), "\n")
    dist1.list[[i]] <- d1
    dist2.list[[i]] <- d2
    pair.list[[i]] <- labels(x[[i]])
    samples.list[[i]] <- labels(d1)
  }
  list(dist1=dist1.list, dist2=dist2.list, pairs=pair.list, 
       samples=samples.list, same.samples=same.samples.list)
}



