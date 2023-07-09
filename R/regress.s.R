# Functions to carry out regressions for s-score regressions
#' @importFrom mgcv gam
#' @importFrom mlr makeRegrTask
#' @importFrom ranger ranger
#' @importFrom tuneRanger tuneRanger
#' @importFrom stats predict lm

regress.s <- function(s_formula, s_learner, boostdf.nuisance.train, boostdf.nuisance.test, boostdf.beta, u_s, U_S, X) {

  scores.df <- cbind(boostdf.nuisance.train[,X,drop=FALSE],u_s=u_s)
  colnames(scores.df)[dim(scores.df)[2]] <- U_S

  if (s_learner == "gam") {
    scores.fit <- gam(s_formula, data=scores.df)
    u_s_hat.train <- scores.fit$fitted
    u_s_hat.test <- predict(scores.fit, boostdf.nuisance.test[,X,drop=FALSE])
    u_s_hat.beta <- predict(scores.fit, boostdf.beta[,X,drop=FALSE])
  } else if (s_learner == "lm") {
    scores.fit <- lm(s_formula, data=scores.df)
    u_s_hat.train <- scores.fit$fitted
    u_s_hat.test <- predict(scores.fit, boostdf.nuisance.test[,X,drop=FALSE])
    u_s_hat.beta <- predict(scores.fit, boostdf.beta[,X,drop=FALSE])
  } else if (s_learner == "randomforest") {
    scores.fit <- ranger(s_formula, scores.df, max.depth=1, num.trees=1) # base learner stumps
    u_s_hat.train <- predict(scores.fit, boostdf.nuisance.train[,X,drop=FALSE])$predictions
    u_s_hat.test <- predict(scores.fit, boostdf.nuisance.test[,X,drop=FALSE])$predictions
    u_s_hat.beta <- predict(scores.fit, boostdf.beta[,X,drop=FALSE])$predictions
  } else {
    stop("Error: Invalid s_learner. s_learner must be gam, lm or randomforest")
  }
  return(list(u_s_hat.train=u_s_hat.train, u_s_hat.test=u_s_hat.test, u_s_hat.beta=u_s_hat.beta))
}
