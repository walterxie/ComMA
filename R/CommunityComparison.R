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
  
  list(proc = proc.list, prot.df=p.result, prot.list=prot.list, pairs=label.matched$pairs, 
       samples=label.matched$samples, same.samples=label.matched$same.samples)
}

#' @details
#' \code{plotProcrustes} plots Procrustes correlations between pairwised communities.  
#' 
#' @param proc.list The list of \code{\link{procrustes}} results.
#' @export
#' @examples 
#' plotProcrustes(procrustes$proc)
#' 
#' @rdname CommunityComparison
plotProcrustes <- function(proc.list, attr.df, colour.id="Elevation") {
  plot.list <- list()
  for(i in 1:length(proc.list)){
    pro <- proc.list[[i]]
    pts <- data.frame(yMDS1 = pro$Yrot[,1], yMDS2 = pro$Yrot[,2], # rotated matrix i.e. d2.mds
                      xMDS1 = pro$X[,1], xMDS2 = pro$X[,2]) # target matrix i.e. d1.mds
    
    pts <- merge(pts, attr.df, by = "row.names")
    
    r1 = acos(pro$rotation[1,1]) # X axis rotation (radians)
    r2 = r1 + (pi/2) # Y axis rotation (radians)
    p <- ggplot(pts) +
      geom_point(aes(x = yMDS1, y = yMDS2, colour = Elevation), shape = 1.5, size = 1.5, alpha = 0.75) + # rotated i.e. d2 (circles)
      geom_point(aes(x = xMDS1, y = xMDS2, colour = Elevation), shape = 2, size = 1, alpha = 0.75) + # target i.e. d1 (triangles)
      scale_shape(solid = FALSE) + xlab("") + ylab("") +
      geom_segment(aes(x = yMDS1, y = yMDS2, xend = xMDS1, yend = xMDS2, colour = Elevation), alpha = 0.75) +
      #arrow = arrow(length = unit(0.2,"cm")), alpha = 0.75) + 
      #geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") + 
      #geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
      #geom_abline(intercept = 0, slope = tan(r1), colour = "grey") + 
      #geom_abline(intercept = 0, slope = tan(r2), colour = "grey") +
      #geom_text(aes(x = xMDS1, y = xMDS2, label = pts.mds$Row.names, colour = Elevation), size = 1.5, vjust = 1.5, alpha = 0.5) + 
      scale_colour_gradientn(colours = c("blue", "orange")) + labs(colour="Elevation (m)") +
      theme(panel.grid = element_blank(), plot.title = element_text(size = 8), 
            plot.margin = unit(c(0.1,0.1,0.1,0), "cm"), legend.key.width = unit(0.65, "cm")) +
      #ggtitle(paste0(letters[[i]], ". ", labels(x[[i]])[[1]], " vs. ", labels(x[[i]])[[2]])) #, ", ", metric, " distance")) #+
      ggtitle(paste0(letters[[i]], ". ", names(proc.list)[[i]]))
    #coord_fixed()
    tmp <- ggplot_gtable(ggplot_build(p))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    
    p <- p + theme(legend.position = "none")
    plot.list[[i]] <- p
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



