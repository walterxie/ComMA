
#' @title
#' Classify sample attribute by taxonomic composition
#' 
#' @details
#' Build a classification model using lasso \code{\link{cv.glmnet}},
#' to predict the certain type of samples (e.g. land types) 
#' given relative abundance of taxonomic compositions (e.g. families) in the samples, 
#' and also to find which taxonomic compositions are important to express this model.
#' 
#' This can be used after 
#' 
#' @param data,attr.data,group.id,percent Refer to \code{\link{enterotype}}.
#' @param alpha,family,nlambda Refer to \code{\link{cv.glmnet}}.
#' @keywords classification
#' @export
#' @examples 
#' cvfit <- classifyByTaxaComp(abun.dist.matrix, attr.data=env, group.id="land.use")
classifyByTaxaComp <- function(data, attr.data, group.id, percent=0.01, 
                               alpha=1, family='multinomial', nlambda = 500, ...) {
  # remove noise
  data <- noise.removal(data, percent=percent)
  x <- t(data)
  y <- attr.data[match(rownames(x), rownames(attr.data)), group.id]
  
  if (nrow(x) != length(y))
    stop("Rows of input matrix x have to be same as response variable y !")
  
  require(glmnet)
  cvfit = cv.glmnet(x, y=as.factor(y), alpha=alpha, family=family, nlambda = nlambda, ...)
  return(cvfit)
}
