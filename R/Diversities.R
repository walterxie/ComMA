# Author: Walter Xie, Alexei Drummond
# Accessed on 9 Sep 2015

######## Jost diversity section #######
# t.community.matrix = t(community.matrix), row is sample

#' @name JostDiversity
#' @title Jost diversity from \pkg{vegetarian} Package
#'
#' @description Data input \strong{t.community.matrix} is 
#' a transposed matrix of community matrix we defined in \pkg{ComMA}.
#' Community matrix from file is a matrix where rows are OTUs or individual species 
#' and columns are sites or samples. See \code{\link{ComMA}}. 
#' 
#' @details \code{diversityTable} creates a data frame 
#' or named vector of Jost diversity measurement.
#' @param t.community.matrix is abundances argument 
#' in \pkg{vegetarian} \code{\link{d}}, 
#' which is a transposed matrix of community matrix, 
#' where rows are plots (Use plots instead of subplots.), 
#' columns are OTUs.
#' @param named.vector Logical, if TRUE, then return 
#' a named vector instead of data.frame.
#' @return 
#' \code{diversityTable} returns a 3x3 data frame: 
#' columns are levels of diversity c("gamma", "alpha", "beta"), 
#' rows are orders of the diversity measure c(0, 1, 2). For example,
#' \tabular{rrrr}{
#'    \tab $q=0$ \tab $q=1$ \tab $q=1$\cr
#'   $D_\\gamma(q)$ \tab 13922.000000 \tab 2501.693162 \tab 601.509610\cr
#'   $D_\\alpha(q)$ \tab 2238.392857 \tab 880.944977 \tab 251.127187\cr
#'   $D_\\beta(q)$ \tab 6.219641 \tab 2.839784 \tab 2.395239 
#' }
#' @export
#' @keywords diversity
#' @examples 
#' diversity.table <- diversityTable(t.community.matrix)
#' 
#' @rdname JostDiversity
diversityTable <- function(t.community.matrix, named.vector=FALSE) { 
  require(vegetarian)
  
  # $ for latex
  diversity <- data.frame(row.names=c("gamma", "alpha", "beta"))
  
  diversity$'$q=0$' <- c(
    d(t.community.matrix,lev="gamma",q=0),
    d(t.community.matrix,lev="alpha",q=0),
    d(t.community.matrix,lev="beta",q=0))
  
  diversity$'$q=1$' <- c(
    d(t.community.matrix,lev="gamma",q=1),
    d(t.community.matrix,lev="alpha",q=1),
    d(t.community.matrix,lev="beta",q=1))
  
  diversity$'$q=2$' <- c(
    d(t.community.matrix,lev="gamma",q=2),
    d(t.community.matrix,lev="alpha",q=2),
    d(t.community.matrix,lev="beta", q=2))
  
  rownames(diversity) <- c("$D_\\gamma(q)$", "$D_\\alpha(q)$", "$D_\\beta(q)$")
  
  if (named.vector) {
    diversity <- unlist(diversity)
    names(diversity) <- c("$D_\\gamma(q=0)$", "$D_\\alpha(q=0)$", "$D_\\beta(q=0)$",
                          "$D_\\gamma(q=1)$", "$D_\\alpha(q=1)$", "$D_\\beta(q=1)$",
                          "$D_\\gamma(q=2)$", "$D_\\alpha(q=2)$", "$D_\\beta(q=2)$")
  }
  return(diversity)
}

#' abundance (reads, gamme0) per sample
#' return 1-column data frame
abundancePerSample <- function(t.community.matrix, hasTotal=TRUE) {
  # gamme0
  perSample <- data.frame(abundance=rowSums(t.community.matrix), stringsAsFactors=FALSE)
  rownames(perSample) <- rownames(t.community.matrix)
  
  if (hasTotal) {
    perSample <- rbind(perSample, sum(perSample[,1]))
    rownames(perSample)[nrow(perSample)] <- "Total"
  }
  
  return(perSample)
}

#' richness (OTUs/species) per sample
#' return 1-column data frame
richnessPerSample <- function(t.community.matrix, hasTotal=TRUE) {
  # richness
  perSample <- data.frame(richness=rowSums(t.community.matrix > 0), stringsAsFactors=FALSE)
  rownames(perSample) <- rownames(t.community.matrix)
  
  if (hasTotal) {
    perSample <- rbind(perSample, sum(perSample[,1]))
    rownames(perSample)[nrow(perSample)] <- "Total"
  }
  
  return(perSample)
}

#' Shannon index (gamma1) per sample
#' return 1-column data frame
shannonPerSample <- function(t.community.matrix, digits = 2) {
  # Shannon
  #  gamma1 <- function(r) d(r,lev="gamma",q=1)
  #  perSample <- data.frame(Shannon=apply(t.community.matrix, 1, gamma1), stringsAsFactors=FALSE)
  perSample <- data.frame(row.names=rownames(t.community.matrix), stringsAsFactors=FALSE)
  require(vegetarian)
  for (i in 1:nrow(t.community.matrix)) 
    perSample[i,1] <- d(t.community.matrix[i,],lev="gamma",q=1)
  
  colnames(perSample)[1] <- "Shannon"
  
  return(round(perSample, digits))
}

#' @details \code{summaryCMPerSample} creates a brief summary table of diversities per sample. 
#' 
#' @param hasTotal If TRUE, then display the Total. Default to TRUE.
#' @param digit The digits of Shannon index. Default to 2.
#' @param sort Sort abundance column by "descending", "ascending", 
#' or do nothing for any other string.
#' @param comma.in.number If TRUE, then add comma to bug numbers. 
#' Default to TRUE.
#' @return 
#' \code{summaryCMPerSample} returns a summary data frame: 
#' columns are types of diversity c("abundance", "richness", "shannon"), 
#' rows are samples. For example,
#' \tabular{rrrr}{
#'    \tab $q=0$ \tab $q=1$ \tab $q=1$\cr
#'   Total \tab 688,733 \tab 16,860 \tab NA\cr
#'   sample1 \tab 32,837 \tab 508 \tab 23.02\cr
#'   sample2 \tab 25,023 \tab 181 \tab 7.9 
#' }
#' @export
#' @keywords diversity
#' @examples 
#' summary.t.cm <- summaryCMPerSample(t.community.matrix)
#' 
#' @rdname JostDiversity
summaryCMPerSample <- function(t.community.matrix, hasTotal=TRUE, digits=2, 
                               sort="descending", comma.in.number=TRUE) {
  # gamme0
  abundance <- abundancePerSample(t.community.matrix, hasTotal) 
  # richness
  richness <- richnessPerSample(t.community.matrix, hasTotal)
  # gamme1
  shannon <- shannonPerSample(t.community.matrix, digits=digits)
  
  if (!all( tolower(rownames(abundance)) == tolower(rownames(richness)) )) 
    stop("data frame rownames do not match !")
  if (!any( tolower(rownames(abundance)[! rownames(abundance) %in% "Total"]) 
            == tolower(rownames(shannon)) )) 
    stop("data frame rownames do not match !")
  
  summary <- data.frame(abundance=abundance, richness=richness, stringsAsFactors=FALSE)
  
  if (hasTotal)
    shannon <- rbind(shannon, c("NA"))
  
  summary<- cbind(summary, shannon)
  if (sort=="descending") {
    summary <- summary[order(summary$abundance, decreasing = TRUE), ]
  } else if (sort=="ascending") {
    summary <- summary[order(summary$abundance), ]
  }# sort=="no"
  if (comma.in.number)
    summary <- format(summary, big.mark=",", scientific=F)
  
  return(summary)
}

