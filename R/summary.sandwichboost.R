#' Summary for a sandwichboost fitted object
#' @description This is a method that prints a \code{sandwich} object fitted by \code{PLR_sandboost()}.
#' @param object a fitted \code{sandwich} object fitted by \code{PLR_sandboost()}.
#' @param ... additional arguments
#' @importFrom stats printCoefmat pt pnorm
#' @export
summary.sandwichboost <- function(object, ...) {
  cat("\nSemiparametric (grouped) partially linear model, with linear component fit by sandwich boosted weighted OLS")
  cat("\nLinear Coefficient:\n")
  printCoefmat(object$coefficients, P.values=TRUE, has.Pvalue=TRUE)
}

#' Print for a sandwichboost fitted object
#' @description This is a method that prints a useful summary of aspects of a \code{sandwich} object fitted by \code{PLR_sandboost()}.
#' @param x a fitted \code{sandwich} object fitted by \code{PLR_sandboost()}.
#' @param ... additional arguments
#' @export
print.sandwichboost <- function(x, ...) {
  cat("\nSemiparametric (grouped) partially linear model, with linear component fit by sandwich boosted weighted OLS")
  cat("\nLinear Coefficient:\n")
  print(x$coefficients, digits=5)
}
