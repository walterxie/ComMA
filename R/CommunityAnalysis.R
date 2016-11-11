
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
#' @param alpha,family,nlambda,... Refer to \code{\link{cv.glmnet}}.
#' @param coef.s Which coeffient to return by \code{\link{coef}}.
#' \code{lambda.min}, as default, is the value of lambda that 
#' gives minimum mean cross-validated error. 
#' \code{lambda.1se} gives the most regularized model such that 
#' error is within one standard error of the minimum.
#' @keywords classification
#' @export
#' @examples 
#' cvfit <- classifyByTaxaComp(abun.dist.matrix, attr.data=env, group.id="land.use")
classifyByTaxaComp <- function(data, attr.data, group.id, percent=0.01, 
                               alpha=1, family='multinomial', nlambda = 500, 
                               coef.s = "lambda.min", return.df=TRUE, ...) {
  # remove noise
  data <- noise.removal(data, percent=percent)
  x <- t(data)
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
