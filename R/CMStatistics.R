# Author: Walter Xie, Andrew Dopheide
# Accessed on 12 Sep 2016

#' @name CMStatistics
#' @title Statistics to one or multipe community matrix
#'
#' @description Statistics to one or multipe community matrix, 
#' such as number of reads, OTUs, etc.  
#' Jost diversity is also included in the summary.
#' 
#' @details \code{summaryCM.Vector} return a named vector of summary of the 
#' community matrix, where \code{community.matrix} can be one column only.
#' The vector is c("reads","OTUs","samples","Shannon","singletons","doubletons").
#' 
#' @param digits The digits to \code{\link{round}} decimal places 
#' if number is not interger. Default to 2.
#' @keywords community matrix
#' @export
#' @examples 
#' summary.cm.vector <- summaryCM.Vector(community.matrix)
#'
#' @rdname CMStatistics
summaryCM.Vector <- function(community.matrix) {
  suppressMessages(require(vegetarian))
  samples <- ncol(community.matrix)
  otus <- nrow(community.matrix)
  reads <- sum(community.matrix)
  if (! is.null(samples)) {
    rs <- rowSums(community.matrix)
    singletons <- sum(rs==1)
    doubletons <- sum(rs==2)
    min.otu.abun <- min(rs)
    max.otu.abun <- max(rs)
    cs <- colSums(community.matrix)
    min.sample.abun <- min(cs)
    max.sample.abun <- max(cs)
  } else {
    samples <- 1
    otus <- length(community.matrix)
    singletons <- NA
    doubletons <- NA
    min.otu.abun <- min(community.matrix)
    max.otu.abun <- max(community.matrix)
    min.sample.abun <- NA
    max.sample.abun <- NA
  }
  shannon <- d(transposeDF(community.matrix),lev="gamma",q=1)
  
  summary.cm <- c(reads, otus, samples, singletons, doubletons, max.otu.abun, 
                  min.otu.abun, max.sample.abun, min.sample.abun)
  names(summary.cm) <- c("reads", "OTUs", "samples", "singletons", "doubletons",
                         "max.OTU.abun","min.OTU.abun","max.sample.abun","min.sample.abun")
  return(summary.cm)
}

#' @details \code{summaryCM} summarizes only single community matrix.
#' 
#' @param has.total If 0, then only return abudence by samples (columns) of community matrix. 
#' If 1, then only return toal abudence. If 2, then return abudence by samples (columns) and total. 
#' Default to 1.
#' @param most.abund The threshold to define the number of the most abundent OTUs.
#' @param pretty.numbers Default to TRUE to make numbers look pretty, 
#' but they will be hard to convert to numeric type.
#' @param x.lab,y.lab,abundance.lab The default text for "sample", "OTU", and "read".
#' @keywords community matrix
#' @export
#' @examples 
#' summary.cm <- summaryCM(community.matrix)
#'
#' @rdname CMStatistics
summaryCM <- function(community.matrix, most.abund, has.total=1, digits=2, pretty.numbers=TRUE,
                      x.lab="sample", y.lab="OTU", abundance.lab="read") {
  summary.row.names <- c(ComMA::getPlural(abundance.lab, y.lab, x.lab),"singletons", "doubletons", 
                       paste("max",y.lab,"abundance",sep="."), paste("min",y.lab,"abundance",sep="."), 
                       paste("max",x.lab,"abundance",sep="."), paste("min",x.lab,"abundance",sep="."))
  summary.cm <- data.frame(row.names = summary.row.names, stringsAsFactors=FALSE, check.names=FALSE)
  if (has.total != 1) {
    for (col.name in colnames(community.matrix)) 
      summary.cm[,col.name] <- ComMA::summaryCM.Vector(community.matrix[,col.name])
  }
  if (has.total > 0) {
    summary.cm[,"total"] <- ComMA::summaryCM.Vector(community.matrix)
    
    if (!missing(most.abund)) {
      if (most.abund > nrow(community.matrix))
        most.abund <- nrow(community.matrix)
      cat("Set most abundent OTUs threshold =", most.abund, ".\n")
      
      community.matrix <- community.matrix[order(rs, decreasing=TRUE),]
      cm <- community.matrix[1:most.abund,]
      col.name <- paste0("most.abund.", most.abund, ".otus)")
      
      summary.cm[,col.name] <- ComMA::summaryCM.Vector(cm)
    }
  }
  if (pretty.numbers)
    summary.cm <- ComMA::prettyNumbers(summary.cm, digits=digits)
  return(summary.cm)
}

#' @details \code{summaryOTUs} returns a data frame of 
#' OTU clustering summary given one or multiple community matrix(matrices).
#' The summary is created by \code{\link{summaryCM.Vector}}.
#' 
#' @param input.list Default to TRUE to unwrap list(...) to 
#' get the actual list if the input is a list of cm. 
#' @param row.names The row names in the summary.
#' Default to "reads", "OTUs", "samples", "singletons", "doubletons",
#' "max.OTU.abun","min.OTU.abun","max.sample.abun","min.sample.abun".
#' @keywords community matrix
#' @export
#' @examples 
#' otu.stats <- summaryOTUs(cm)
#'
#' @rdname CMStatistics
summaryOTUs <- function(..., digits=2, input.list=FALSE, pretty.numbers=TRUE,
                        x.lab="sample", y.lab="OTU", abundance.lab="read") {
  cm.list <- list(...)
  # if input a list of cm, then unwrap list(...) to get the actual list
  if (input.list && typeof(cm.list[[1]]) == "list")
    cm.list <- cm.list[[1]]
  # validation
  if (! (typeof(cm.list) == "list" && length(cm.list) > 0) )
    stop("Invaild input: at least one community matrix is required ! length(...) = ", length(cm.list))
  
  cat("Summarize OTUs on", length(cm.list), "data sets.\n") 
  
  summary.row.names <- c(ComMA::getPlural(abundance.lab, y.lab, x.lab),"singletons", "doubletons", 
                         paste("max",y.lab,"abundance",sep="."), paste("min",y.lab,"abundance",sep="."), 
                         paste("max",x.lab,"abundance",sep="."), paste("min",x.lab,"abundance",sep="."))
  otu.stats <- data.frame(row.names = summary.row.names, stringsAsFactors=FALSE, check.names=FALSE)
  for (data.id in 1:length(cm.list)) {
    col.name <- names(cm.list)[data.id]
    if (is.null(col.name))
      col.name <- paste("data",data.id,sep=".")
    # otus.row.names
    otu.stats[,col.name] <- ComMA::summaryCM.Vector(cm.list[[data.id]])
  }
  otu.stats["non.singletons",] <- as.numeric(otu.stats["OTUs",])-as.numeric(otu.stats["singletons",])
  
  if (pretty.numbers) 
    otu.stats <- ComMA::prettyNumbers(otu.stats, digits=digits)
  return(otu.stats)
}

#' @details \code{summaryDiversity} returns a data frame of 
#' OTU clustering summary given a list of community matrix.
#' The community matrix is transposed to an input of 
#' \code{\link{d}} {vegetarian}, and the summary is 
#' created by \code{\link{diversityTable}}.
#' 
#' @param row.names The row names in the summary.
#' Default to "$^0D_\\gamma$","$^0D_\\alpha$","$^0D_\\beta$",
#' "$^1D_\\gamma$","$^1D_\\alpha$", "$^1D_\\beta$",
#' "$^2D_\\gamma$","$^2D_\\alpha$","$^2D_\\beta$".
#' @param row.order The same row indices of \code{row names}, 
#' but by a given order. Default to c().
#' @keywords community matrix
#' @export
#' @examples 
#' div.stats <- summaryDiversity(cm, row.order=c(2,5,8,3,6,9,1,4,7))
#'
#' @rdname CMStatistics
summaryDiversity <- function(..., row.order=c(), digits=2, input.list=FALSE, pretty.numbers=TRUE, verbose=TRUE, 
        row.names=c("$^0D_\\gamma$","$^0D_\\alpha$","$^0D_\\beta$","$^1D_\\gamma$","$^1D_\\alpha$",
                    "$^1D_\\beta$","$^2D_\\gamma$","$^2D_\\alpha$","$^2D_\\beta$")) {
  cm.list <- list(...)
  # if input a list of cm, then unwrap list(...) to get the actual list
  if (input.list && typeof(cm.list[[1]]) == "list")
    cm.list <- cm.list[[1]]
  # validation
  if (! (typeof(cm.list) == "list" && length(cm.list) > 0) )
    stop("Invaild input: at least one community matrix is required ! length(...) = ", length(cm.list))
  
  cat("Summarize Jost diversities on", length(cm.list), "data sets.\n") 
  
  div.stats <- data.frame(stringsAsFactors=FALSE, check.names=FALSE)
  for (data.id in 1:length(cm.list)) {
    t.cm <- ComMA::transposeDF(cm.list[[data.id]])
    # gamma(q=0), alpha(q=0), beta(q=0), gamma(q=1), ..., beta(q=2)
    diversity.v <- ComMA::diversityTable(t.cm, named.vector=T)
    i.vec <- 1:length(diversity.v)
    if (length(row.order) == length(row.names)) {
      if (! all( 1:length(row.names) == sort(row.order) ) )
        stop("Invalid vector row.order !")
      # reorder to match div.row.names
      i.vec <- row.order
    } 
    
    col.name <- names(cm.list)[data.id]
    if (is.null(col.name))
      col.name <- paste("data",data.id,sep=".")

    for (i in i.vec) 
      div.stats[row.names[i], col.name] <- diversity.v[i]
    
    if (verbose)
      cat("Add column", col.name, "to data set", data.id, ".\n")
  }
  
  if (pretty.numbers) 
    div.stats <- ComMA::prettyNumbers(div.stats, digits=digits)
  return(div.stats)
}


#' @details \code{summaryTaxaGroup} returns two data frames of 
#' taxonomic composition of given one or multiple merged 
#' community matrix(matrices) with taxa table, 
#' which is created by \code{\link{mergeCMTaxa}}.
#' The 1st data frame \code{otus} summarizes the OTUs 
#' assigned to each taxa group.
#' The 2nd data frame \code{rank.count} summarizes the \code{count.rank} 
#' assigned to each taxa group.
#' 
#' @param unclassified Refere to \code{\link{assignTaxaByRank}}, 
#' default to 3 to remove every rows containing "unclassified". 
#' @param taxa.group The row names in the summary.
#' Default to "ARCHAEA", "BACTERIA", "CHROMISTA", "PROTOZOA", 
#' "FUNGI", "PLANTAE", "ANIMALIA", "EUKARYOTA".
#' @param group.rank The rank of given taxa groups.
#' @param count.rank The lower rank to be counted for each taxa group. 
#' Set NA to ignore this count.
#' @keywords community matrix
#' @export
#' @examples 
#' ta.gr.stats <- summaryTaxaGroup(cm.taxa)
#' ta.gr.stats$otus
#' ta.gr.stats$rank.count 
#'
#' @rdname CMStatistics
summaryTaxaGroup <- function(..., input.list=FALSE, unclassified=3, pretty.numbers=TRUE,
          taxa.group=c("ARCHAEA", "BACTERIA", "CHROMISTA", "PROTOZOA", "FUNGI", "PLANTAE", "ANIMALIA"), 
          group.rank="kingdom", count.rank="phylum") {
  cm.taxa.list <- list(...)
  # if input a list of cm, then unwrap list(...) to get the actual list
  if (input.list && typeof(cm.taxa.list[[1]]) == "list")
    cm.taxa.list <- cm.taxa.list[[1]]
  # validation
  if (! (typeof(cm.taxa.list) == "list" && length(cm.taxa.list) > 0) )
    stop("Invaild input: at least one community matrix is required ! length(...) = ", length(cm.taxa.list))
  
  cat("Summarize OTUs by taxa group on", length(cm.taxa.list), "data sets.\n") 
  
  # data frame for statistics
  tg.reads <- data.frame(row.names = taxa.group, stringsAsFactors=FALSE, check.names=FALSE)
  tg.otus <- data.frame(row.names = taxa.group, stringsAsFactors=FALSE, check.names=FALSE)
  tg.rank.count <- NA
  if (!is.na(count.rank))
    tg.rank.count <- data.frame(row.names = taxa.group, stringsAsFactors=FALSE, check.names=FALSE)
  count.rank.df.list <- list()
  
  for (data.id in 1:length(cm.taxa.list)) {
    col.name <- names(cm.taxa.list)[data.id]
    if (is.null(col.name))
      col.name <- paste("data",data.id,sep=".")
    
    for (taxa.id in 1:length(taxa.group)) {
      cat("Data set", col.name, "taxonomic group", taxa.group[taxa.id], ".\n") 
      
      cm.taxa.sub <- ComMA::subsetTaxaTable(cm.taxa.list[[data.id]], taxa.group=taxa.group[taxa.id], rank=group.rank)
      
      if (nrow(cm.taxa.sub) < 1) {
        tg.reads[taxa.group[taxa.id],col.name] <- 0
        tg.otus[taxa.group[taxa.id],col.name] <- 0
        if (!is.na(count.rank))
          tg.rank.count[taxa.group[taxa.id],col.name] <- 0
      } else {
        taxa.assign.otu <- ComMA::assignTaxaByRank(cm.taxa.sub, unclassified=unclassified, aggre.FUN=function(x) sum(x>0))
        taxa.assign.reads <- ComMA::assignTaxaByRank(cm.taxa.sub, unclassified=unclassified)
        
        tg.otus[taxa.group[taxa.id],col.name] <- sum(taxa.assign.otu[[group.rank]])
        tg.reads[taxa.group[taxa.id],col.name] <- sum(taxa.assign.reads[[group.rank]])
        if (!is.na(count.rank)) {
          tg.rank.count[taxa.group[taxa.id],col.name] <- nrow(taxa.assign.otu[[count.rank]])
          # record taxa in count.rank
          df.name <- paste(col.name, taxa.group[taxa.id], sep = ".")
          count.rank.df.list[[df.name]] <- taxa.assign.otu[[count.rank]]
        }
      }
    }
  }
  
  if (pretty.numbers) {
    tg.otus <- ComMA::prettyNumbers(tg.otus, digits = 0)
    tg.reads <- ComMA::prettyNumbers(tg.reads, digits = 0)
    if (!is.na(count.rank))
      tg.rank.count <- ComMA::prettyNumbers(tg.rank.count, digits = 0)
  }
  
  list(otus=tg.otus, rank.count=tg.rank.count, reads=tg.reads, taxa.group=taxa.group, group.rank=group.rank,
      count.rank=count.rank, count.rank.df.list=count.rank.df.list)
}

