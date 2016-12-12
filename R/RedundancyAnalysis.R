# A Caution Regarding Rules of Thumb for Variance Inflation Factors
# http://link.springer.com/article/10.1007/s11135-006-9018-6

#' @name RDA
#' @title Redundancy Analysis
#'
#' @description Redundancy Analysis using \pkg{vegan} 
#' \code{\link{capscale}}.
#' 
#' @details \code{proceedRDA} makes Constrained Analysis 
#' of Principal Coordinates for 
#' eDNA data sets given environmental variables.
#' 
#' @param tcm.or.dist A transposed community matrix 
#' or dist object of distances between samples. 
#' Rows are samples.
#' @param env The enviornmental meta-data, where rows are samples, 
#' and they must be same as rownames(tcm.or.dist) inlcuding order.
#' In addition, make sure rownames (enviornmental variables) 
#' are valid to \code{\link{formula}}.
#' @param verbose More details. Default to TRUE.
#' @keywords rda
#' @export
#' @examples 
#' rda <- proceedRDA(tcm.or.dist, env)
#' 
#' @rdname RDA
proceedRDA <- function(tcm.or.dist, env, verbose=TRUE) {
  if (! all( tolower(rownames(env)) == tolower(rownames(as.matrix(tcm.or.dist))) ) ) 
    stop("Site names in community matrix and environmental data file not matched !")
  
  require(vegan)
  #require(VIF)
  # Constrained ordination ------------------------------------------------------
  rda_table <- data.frame(row.names=c("Constrained","Unconstrained"))
  anova_table <- data.frame(row.names=colnames(env))
  
  # Distance-based redundancy analysis, using capscale
  # Constrained Analysis of Principal Coordinates (CAP) is an ordination method similar to Redundancy Analysis (rda).
  # DB-RDA, empty model
  rda_0 <- capscale(tcm.or.dist ~ 1, env, distance = "jaccard")
  
  # DB-RDA, maximal model (bad idea - only use for auto model building)
  rda_1 <- capscale(tcm.or.dist ~ ., env, distance = "jaccard")
  if (verbose) head(summary(rda_1))
  # sp = species scores, wa = site scores, bp = biplot arrows, lc = linear constraints 
  #	plot(rda_1, display = c("wa", "bp")) # Note correlation of biplot arrows
  
  # Variance inflation factor - indicates highly correlated variables
  if (verbose) print(vif.cca(rda_1))
  
  rda_table$Inertia <- c(round(rda_1$CCA$tot.chi, 3), round(rda_1$CA$tot.chi, 3))
  rda_table$Proportion <- c(rda_1$CCA$tot.chi/rda_1$tot.chi, rda_1$CA$tot.chi/rda_1$tot.chi)
  rda_table$Proportion <- percent(rda_table$Proportion) # %
  
  constrained_inertia <- c()
  constrained_proportion <- c()
  # Test each variable individually
  for (i in 1:length(colnames(env))) {
    rda_individual <- capscale(formula = as.formula(paste("tcm.or.dist", colnames(env)[i], sep=" ~ ")), 
                               env, distance = "jaccard")
    constrained_inertia <- c(constrained_inertia, rda_individual$CCA$tot.chi)
    constrained_proportion <- c(constrained_proportion, rda_individual$CCA$tot.chi/rda_individual$tot.chi)
  }
  
  anova_table$Inertia <- constrained_inertia
  anova_table$Proportion <- percent(constrained_proportion) # %
  
  #Compute all the single terms in the scope argument that can be added to or dropped from the model, 
  #fit those models and compute a table of the changes in fit.
  add_1 <- add1(rda_0, scope=formula(rda_1), test="perm")
  anova_table$Pr <- add_1$Pr[-1] # 1st row is <none>
  colnames(anova_table)[3] <- "Pr($>$F)"
  
  # Build model after stepwise removal of collinear variables (vif >= 10; requires vif_function.R) 
  # variance inflation factor (VIF) quantifies the severity of multicollinearity in an ordinary least squares regression analysis. 
  env_reduced <- vif_func(in_frame = env)
  if (verbose) print(env_reduced) # Remaining variables
  
  # Build model automatically from reduced variable set
  # (Unsure how to pass env_reduced variables to capscale formula; paste() doesn't work...)
  #	rda_reduced <- capscale(tcm.or.dist ~ slope.degree + Mean.Temp + Northness + Eastness + 
  #							pH + C.N.ratio + NO3.N + NH4.N + Olsen.P, env, distance = "jaccard")
  rda_reduced <- capscale(formula = as.formula(paste("tcm.or.dist", paste(env_reduced, collapse=" + "), sep=" ~ ")), 
                          env, distance = "jaccard")
  if (verbose) head(summary(rda_reduced))
  anova_reduced <- anova(rda_reduced, by = "terms")
  
  rda_table$Inertia.R <- c(round(rda_reduced$CCA$tot.chi, 3), round(rda_reduced$CA$tot.chi, 3))
  rda_table$Proportion.R <- c(rda_reduced$CCA$tot.chi/rda_reduced$tot.chi, rda_reduced$CA$tot.chi/rda_reduced$tot.chi)
  rda_table$Proportion.R <- percent(rda_table$Proportion.R) # %
  
  anova_table$Reduced <- is.element(rownames(anova_table), rownames(anova_reduced))
  anova_table$Reduced[which(anova_table$Reduced==T)] <- anova_reduced$Pr[-length(anova_reduced$Pr)]
  colnames(anova_table)[4] <- "Reduced Pr($>$F)"
  
  # Choose a Model by Permutation Tests in Constrained Ordination using forward model selection
  rda_reduced_f <- ordistep(rda_0, scope = formula(rda_reduced), direction = "forward", permutations = 3999)
  
  rda_forward <- capscale(formula = as.formula(rda_reduced_f$call), data = env, distance = "jaccard")
  if (verbose) head(summary(rda_forward))
  anova_forward <- anova(rda_forward, by = "terms")
  
  rda_table$Inertia.F <- c(round(rda_forward$CCA$tot.chi, 3), round(rda_forward$CA$tot.chi, 3))
  rda_table$Proportion.F <- c(rda_forward$CCA$tot.chi/rda_forward$tot.chi, rda_forward$CA$tot.chi/rda_forward$tot.chi)
  rda_table$Proportion.F <- percent(rda_table$Proportion.F) # %
  
  anova_table$Forward <- is.element(rownames(anova_table), rownames(anova_forward))
  anova_table$Forward[which(anova_table$Forward==T)] <- anova_forward$Pr[-length(anova_forward$Pr)]
  colnames(anova_table)[5] <- "Forward Pr($>$F)"
  
  # Choose a Model by Permutation Tests in Constrained Ordination using backward model selection
  rda_reduced_b <- ordistep(rda_reduced, scope = formula(rda_0), direction = "backward", permutations = 3999)
  
  rda_backward <- capscale(formula = as.formula(rda_reduced_b$call), data = env, distance = "jaccard")
  if (verbose) head(summary(rda_backward))
  anova_backward <- anova(rda_backward, by = "terms")
  
  rda_table$Inertia.B <- c(round(rda_backward$CCA$tot.chi, 3), round(rda_backward$CA$tot.chi, 3))
  rda_table$Proportion.B <- c(rda_backward$CCA$tot.chi/rda_backward$tot.chi, rda_backward$CA$tot.chi/rda_backward$tot.chi)
  rda_table$Proportion.B <- percent(rda_table$Proportion.B) # %
  
  anova_table$Backward <- is.element(rownames(anova_table), rownames(anova_backward))
  anova_table$Backward[which(anova_table$Backward==T)] <- anova_backward$Pr[-length(anova_backward$Pr)]
  colnames(anova_table)[6] <- "Backward Pr($>$F)"
  
  anova_table[anova_table == 0] <- ""
  anova_table$Proportion <- gsub("%", "\\\\%", anova_table$Proportion)
  
  rda_table$Proportion <- gsub("%", "\\\\%", rda_table$Proportion)
  rda_table$Proportion.R <- gsub("%", "\\\\%", rda_table$Proportion.R)
  rda_table$Proportion.F <- gsub("%", "\\\\%", rda_table$Proportion.F)
  rda_table$Proportion.B <- gsub("%", "\\\\%", rda_table$Proportion.B)
  
  # Return a list 
  list( reduced = rda_reduced, forward = rda_forward, backward = rda_backward,
        anova.summary=anova_table, model.summary=rda_table)
}

#' @details \code{preprocessRDA} subsets variables and samples to 
#' be included in analysis \code{proceedRDA}.
#' It is required that samples in community matrix must match 
#' the enviornmental meta-data, including the order. 
#' Teh function does not handle dist object of distances.
#' 
#' @param cm A community matrix. 
#' @param is.transposed If TRUE, then the community matrix is already
#' transposed to be the valid input of \code{\link{vegdist}}.  
#' Default to FASLE.
#' @param rm.samples Remove specified samples in a vector, 
#' it can be a keyword shared in sample names.
#' The vector will convert to a string separated by '|' to multi-samples. 
#' Default to empty vector to do nothing.
#' @param min.abund Exclude any samples with excessively 
#' low abundance. Defaul 0 to exaclude none.
#' The code is \code{tcm <- tcm[rowSums(tcm) > min.abund, ]}. 
#' @param sel.env.var The vector of selected environmental variables, 
#' which can be colnames(env) or their indices. 
#' Defaul to an empty vector to choose all variables.
#' @param log.var,log.base It normally needs log transform to soil chemistry variables.
#' Use \code{\link{plotCorrelations}} to visualize variables and determine 
#' whether log transform should be applied. Default to no log transform. 
#' @export
#' @examples 
#' tcm.env <- preprocessRDA(tcm, env, is.transposed=T)
#' # Note colSums(cm) are based on samples
#' tcm.env <- preprocessRDA(cm, env, rm.samples=c("CM30b51","CM30b58"), min.abund=mean(colSums(cm))*0.025, sel.env.var=c(4,5,8,9,14:22))
#' 
#' @rdname RDA
preprocessRDA <- function(cm, env, is.transposed=FALSE, rm.samples=c(), min.abund=0, 
                          sel.env.var=c(), log.var=c(), log.base=2) {
  if (!is.transposed) {
    print(ComMA::summaryCM(cm))
    tcm <- ComMA::transposeDF(cm)
  } else {
    print(ComMA::summaryCM(t(cm)))
    tcm <- cm
  } 
  
  # remove specified samples, it can be keywords.
  if (length(rm.samples) > 0) {
    rm <- paste(rm.samples, collapse = "|")
    env <- env[!grepl(rm, rownames(env), ignore.case = T), ]
    n.samples <- nrow(tcm)
    tcm <- tcm[!grepl(rm, rownames(tcm), ignore.case = T), ]
    cat("Drop", n.samples-nrow(tcm), "samples containing : ", paste(rm.samples, collapse = ","), "\n") 
  }
  # exclude any samples with excessively low abundance
  if (min.abund > 0) {
    print(summary(rowSums(tcm)))
    n.samples <- nrow(tcm)
    tcm <- tcm[rowSums(tcm) > min.abund, ]
    tcm <- tcm[, colSums(tcm) > 0] # Exclude any empty col 
    cat("Drop", n.samples-nrow(tcm), "samples with low abundance <= ", min.abund, "\n") 
  }
  
  # select environmental variables
  if (length(sel.env.var) > 0) {
    env <- env[, sel.env.var]
    cat("Select", length(sel.env.var), "environmental variables : ", paste(colnames(env), collapse = ","), "\n") 
  }
  # select environmental variables
  if (length(log.var) > 0) {
    env[,log.var] <- log(env[,log.var], log.base)
    env[env == "-Inf"] <- 0 # Replace inf with zero
    cat("Log transform", length(log.var), "variables : ", paste(colnames(env[,log.var]), collapse = ","), 
        "at base", log.base, "\n") 
  }
  
  if (! all( tolower(rownames(env)) == tolower(rownames(tcm)) ) ) { # Rows don't match
    both <- intersect(rownames(tcm), rownames(env)) # Matching rows
    if (length(both)<1) 
      stop("Invalid transposed community matrix and enviornmental meta-data : no row names matching !")
    # Subset to matching rows
    if (length(both) < nrow(tcm)) 
      cat("Drop samples not present in tcm : ", paste(setdiff(rownames(tcm), both), collapse = ","), "\n")
    tcm <- tcm[match(both, rownames(tcm)),]
    if (length(both) < nrow(env)) 
      cat("Drop samples not present in env : ", paste(setdiff(rownames(env), both), collapse = ","), "\n")
    env <- env[match(both, rownames(env)),] 
  }
  cat("After preprocessing, tcm has", nrow(tcm), "samples,", ncol(tcm), "OTUs; env has", 
      nrow(env), "samples,", ncol(env), "enviornmental variables.\n\n")
  # replace invalid char for formula  
  colnames(env) <- gsub("-", ".", colnames(env))
  colnames(env) <- gsub("/", ".", colnames(env))
  list(tcm=tcm, env=env)
}

#' @details \code{plotCorrelations} plots numeric variables (columns).
#' 
#' Tip: use \code{\link{"\%<a-\%"}} in \pkg{pryr} to save plots.
#' 
#' @param df.numeric The data frame or matrix containing 
#' numeric variables (columns) to plot.
#' @param corr.gram Logical, if use \code{\link{corrgram}} 
#' instead of \code{\link{plot}}. 
#' @export
#' @examples 
#' # before RDA
#' require(pryr)
#' p %<a-% plotCorrelations(tcm.env$env)
#' 
#' @rdname RDA
plotCorrelations <- function(df.numeric, corr.gram=FALSE, cex.axis = 0.75, 
                             cex.cor = 0.9, col = "#333333") {
  if (corr.gram) 
    corrgram(df.numeric, gap = 0, lower.panel = panel.pts, upper.panel = panel.conf, 
             cex.axis = cex.axis, cex.cor = cex.cor, col = col)
  else 
    plot(df.numeric, gap = 0, lower.panel = panel.smooth, upper.panel = panel.conf, 
         cex.axis = cex.axis, cex.cor = cex.cor, col.smooth = "purple", col = col)
}

#' @details \code{printXTable.RDA} prints \code{\link{xtable}} given rda results.
#' 
#' @param rda The list of results from \code{proceedRDA}.
#' @param matrix.name The string to locate the matrix from its file name. 
#' Only used for table name and label here.
#' @param taxa.group The taxonomic group. Only used for table name and label here. 
#' @param table.file If NULL, then print the results to console, 
#' otherwise print them to the file. Default to NULL.
#' @examples 
#' printXTable.RDA(rda, table.file=NULL, matrix.name="16S", taxa.group="BACTERIA")
#' 
#' @rdname RDA
printXTable.RDA <- function(rda, table.file=NULL, invalid.char=FALSE, matrix.name="", taxa.group="") {
  
  ComMA::printXTable(rda$anova.summary, invalid.char=invalid.char,
              caption = paste("Distance-based redundancy analysis and their ANOVA tests 
                              in each step for the eDNA biodiversity data sets", matrix.name, taxa.group), 
              label = paste("tab:rdaAnova", matrix.name, taxa.group, sep = ":"), file=table.file)
  
  ComMA::printXTable(rda$model.summary, invalid.char=invalid.char,
              caption = paste("The constrained and unconstrained inertia changes during 
          distance-based redundancy analysis for the eDNA biodiversity data sets", matrix.name, taxa.group), 
              label = paste("tab:rda", matrix.name, taxa.group, sep = ":"), file=table.file)
}
