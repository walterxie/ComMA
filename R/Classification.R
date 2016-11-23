
#' @name classification
#' @title Classify samples by taxonomic composition
#'
#' @description
#' Build a classification model to classify or predict 
#' the certain type of samples (e.g. land types).
#' 
#' @details
#' Create a conditional inference tree \code{\link{ctree}} for the classification of samples.
#' Ctree uses a significance test procedure in order to select variables 
#' instead of selecting the variable that maximizes an information measure (e.g. Gini coefficient).
#' \url{http://stats.stackexchange.com/questions/12140/conditional-inference-trees-vs-traditional-decision-trees}.
#' 
#' @param comm A community matrix, which can be either abundance or relative abundance.
#' @param levels Levels to order data by \code{group.id}.
#' @keywords classification
#' @export
#' @examples 
#' model <- ctreeClassification(comm, attr.data=env, group.id="land.use")
#' plot(model$ctree)
#' 
#' @rdname classification
ctreeClassification <- function(comm, attr.data, group.id, levels=c(), ...) {
  require(party)
  input.dat <- as.data.frame(t(comm))
  input.dat[,group.id] <- attr.data[match(rownames(input.dat), rownames(attr.data)), group.id]
  if (length(levels) < 1)
    levels <- sort(unique(input.dat[,group.id]))
  input.dat[,group.id] <- factor(input.dat[,group.id], levels=levels)
  
  cat("Input community has", nrow(input.dat), "samples, and the target", 
      group.id, "has", length(unique(input.dat[,group.id])), "categories.")
  
  output.tree <- ctree(as.formula(paste(group.id, ".", sep=" ~ ")), data = input.dat, ...)
  return(ctree=output.tree, data=input.dat, levels=levels)
}

#' @details
#' Build a classification model using lasso \code{\link{cv.glmnet}},
#' given relative abundance of taxonomic compositions (e.g. families) in the samples, 
#' and also to find which taxonomic compositions are important.
#' 
#' @param relative.abund,attr.data,percent Refer to \code{\link{enterotype}}.
#' @param group.id The column name in \code{attr.data} contains the known groups to 
#' compare with enterotypes.
#' @param alpha,family,nlambda,... Refer to \code{\link{cv.glmnet}}.
#' @param coef.s Which coeffient to return by \code{\link{coef}}.
#' \code{lambda.min}, as default, is the value of lambda that 
#' gives minimum mean cross-validated error. 
#' \code{lambda.1se} gives the most regularized model such that 
#' error is within one standard error of the minimum.
#' @keywords classification
#' @export
#' @examples 
#' cvfit <- lassoClassification(relative.abund, attr.data=env, group.id="land.use")
#' 
#' @rdname classification
lassoClassification <- function(relative.abund, attr.data, group.id, percent=0.01, 
                                alpha=1, family='multinomial', nlambda = 500, 
                                coef.s = "lambda.min", return.df=TRUE, ...) {
  # remove noise
  relative.abund <- noise.removal(relative.abund, percent=percent)
  x <- t(relative.abund)
  y <- attr.data[match(rownames(x), rownames(attr.data)), group.id]
  
  if (nrow(x) != length(y))
    stop("Rows of input matrix x have to be same as response variable y !")
  
  require(glmnet)
  cvfit = cv.glmnet(x, y=as.factor(y), alpha=alpha, family=family, nlambda = nlambda, ...)
  # 
  coef.list <- coef(cvfit, s = coef.s)
  model <- do.call("cbind", coef.list)
  colnames(model) <- names(coef.list)
  if (return.df) {
    model <- as.data.frame(as.matrix(model))
    model <- model[-1,]                  # rm 1st row "(Intercept)" 
    model <- model[rowSums(model) != 0,]  # rm empty rows
  }
  
  list(fit=cvfit, model=model, return.df=return.df)
}

