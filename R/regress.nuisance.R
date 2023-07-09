# Functions to carry out regressions for nuisance functions.
#' @importFrom mgcv gam
#' @importFrom glmnet cv.glmnet
#' @importFrom mlr makeRegrTask
#' @importFrom ranger ranger
#' @importFrom tuneRanger tuneRanger
#' @importFrom stats lm predict sigma

# Regression for m(x)=E[D|X=x].
regress.m <- function(m_formula, m_learner, data.nuisance, data.beta) {
  D <- all.vars(m_formula)[1]
  X <- all.vars(m_formula)[-1]
  if (m_learner == "gam") {
    nuisance_m_fit <- gam(m_formula, data=data.nuisance)
    m_residuals.nuisance <- data.nuisance[,D] - nuisance_m_fit$fitted
    m_residuals.beta <- data.beta[,D] - predict(nuisance_m_fit, data.beta[,X,drop=FALSE])
  } else if (m_learner == "lm") {
    nuisance_m_fit <- lm(m_formula, data.nuisance)
    m_residuals.nuisance <- data.nuisance[,D] - nuisance_m_fit$fitted
    m_residuals.beta <- data.beta[,D] - predict(nuisance_m_fit, data.beta[,X,drop=FALSE])
  } else if (m_learner == "randomforest") {
    stuff <- data.nuisance[, (names(data.nuisance) %in% c(D, X))]
    m.task <- makeRegrTask(data=stuff, target="D")
    res <- tuneRanger(m.task, iters = 70, iters.warmup = 30, time.budget = NULL, num.threads = NULL, num.trees = 500, parameters = list(replace = FALSE, respect.unordered.factors="order"), tune.parameters = c("min.node.size"), show.info=FALSE)
    nuisance_m_fit <- ranger(m_formula, data.nuisance, min.node.size=res$recommended.pars$min.node.size, num.trees=500)
    m_residuals.nuisance <- data.nuisance[,D] - predict(nuisance_m_fit, data.nuisance[,X,drop=FALSE])$predictions
    m_residuals.beta <- data.beta[,D] - predict(nuisance_m_fit, data.beta[,X,drop=FALSE])$predictions
  } else if (m_learner == "lasso") {
    nuisance_m_fit <- cv.glmnet(as.matrix(data.nuisance[,X]), as.matrix(data.nuisance[,D]))
    m_residuals.nuisance <- data.nuisance[,D] - predict(nuisance_m_fit, newx = as.matrix(data.nuisance[,X]), s = "lambda.min")
    colnames(m_residuals.nuisance) <- NULL
    m_residuals.beta <- data.beta[,D] - predict(nuisance_m_fit, newx = as.matrix(data.beta[,X]), s = "lambda.min")
    colnames(m_residuals.beta) <- NULL
  } else if (m_learner == "ridge") {
    nuisance_m_fit <- cv.glmnet(as.matrix(data.nuisance[,X]), as.matrix(data.nuisance[,D]), alpha=0)
    m_residuals.nuisance <- data.nuisance[,D] - predict(nuisance_m_fit, newx = as.matrix(data.nuisance[,X]), s = "lambda.min")
    colnames(m_residuals.nuisance) <- NULL
    m_residuals.beta <- data.beta[,D] - predict(nuisance_m_fit, newx = as.matrix(data.beta[,X]), s = "lambda.min")
    colnames(m_residuals.beta) <- NULL
  } else {
    stop("Error: Invalid m_learner. m_learner must be gam, randomforest, lasso, ridge, or lm.")
  }
  return(list(m_residuals.nuisance=m_residuals.nuisance, m_residuals.beta=m_residuals.beta))
}

# Regression for l(x)=E[Y|X=x].
regress.l <- function(l_formula, l_learner, data.nuisance, data.beta) {
  Y <- all.vars(l_formula)[1]; X <- all.vars(l_formula)[-1]
  if (l_learner == "gam") {
    nuisance_l_fit <- gam(l_formula, data=data.nuisance)
    l_residuals.nuisance <- data.nuisance[,Y] - nuisance_l_fit$fitted
    l_residuals.beta <- data.beta[,Y] - predict(nuisance_l_fit, data.beta[,X,drop=FALSE])
  } else if (l_learner == "lm") {
    nuisance_l_fit <- lm(l_formula, data.nuisance)
    l_residuals.nuisance <- data.nuisance[,Y] - nuisance_l_fit$fitted
    l_residuals.beta <- data.beta[,Y] - predict(nuisance_l_fit, data.beta[,X,drop=FALSE])
  } else if (l_learner == "randomforest") {
    stuff <- data.nuisance[, names(data.nuisance) %in% c(Y,X)]
    l.task <- makeRegrTask(data=stuff, target=Y)
    res <- tuneRanger(l.task, iters = 70, iters.warmup = 30, time.budget = NULL, num.threads = NULL, num.trees = 500, parameters = list(replace = FALSE, respect.unordered.factors="order"), tune.parameters = c("min.node.size"), show.info=FALSE)
    nuisance_l_fit <- ranger(l_formula, data.nuisance, min.node.size=res$recommended.pars$min.node.size, num.trees=500)
    l_residuals.nuisance <- data.nuisance[,Y] - predict(nuisance_l_fit, data.nuisance[,X,drop=FALSE])$predictions
    l_residuals.beta <- data.beta[,Y] - predict(nuisance_l_fit, data.beta[,X,drop=FALSE])$predictions
  } else if (l_learner == "lasso") {
    nuisance_l_fit <- cv.glmnet(as.matrix(data.nuisance[,X]), as.matrix(data.nuisance[,Y]))
    l_residuals.nuisance <- data.nuisance[,Y] - predict(nuisance_l_fit, newx = as.matrix(data.nuisance[,X]), s = "lambda.min")
    colnames(l_residuals.nuisance) <- NULL
    l_residuals.beta <- data.beta[,Y] - predict(nuisance_l_fit, newx = as.matrix(data.beta[,X]), s = "lambda.min")
    colnames(l_residuals.beta) <- NULL
  } else if (l_learner == "ridge") {
    nuisance_l_fit <- cv.glmnet(as.matrix(data.nuisance[,X]), as.matrix(data.nuisance[,Y]), alpha=0)
    l_residuals.nuisance <- data.nuisance[,Y] - predict(nuisance_l_fit, newx = as.matrix(data.nuisance[,X]), s = "lambda.min")
    colnames(l_residuals.nuisance) <- NULL
    l_residuals.beta <- data.beta[,Y] - predict(nuisance_l_fit, newx = as.matrix(data.beta[,X]), s = "lambda.min")
    colnames(l_residuals.beta) <- NULL
  } else {
    stop("Error: Invalid l_learner. l_learner must be gam, randomforest, lasso, ridge, or lm.")
  }
  return(list(l_residuals.nuisance=l_residuals.nuisance, l_residuals.beta=l_residuals.beta))
}
