#' Cross Validation of the sandwich loss for selection of number of boosting iterations
#'
#' @param l_formula a two-sided formula object describing the regression model for \eqn{E[Y|X=x]}.
#' @param l_learner a string specifying the regression method to fit the regression of \eqn{Y} on \eqn{X} as given by \code{l_formula}.
#' @param m_formula a two-sided formula object describing the regression model for \eqn{E[D|X=x]}.
#' @param m_learner a string specifying the regression method to fit the regression of \eqn{D} on \eqn{X} as given by \code{m_formula}.
#' @param s_formula a two-sided formula object describing the regression model for a single boosting iteration for \eqn{s}-function sandwich boosting.
#' @param s_learner a string specifying the regression method to fit the a single boosting iteration for \eqn{s}-function sandwich boosting as given by \code{l_formula}.
#' @param proxyCCF a string specifying which working correlation parametrisation to use. Must take the form "\code{equicorr}", "\code{autoreg}" or "\code{nested}".
#' @param data a data frame containing the variables for the grouped PLR model. Note that group ID must be given by the column \code{id}.
#' @param K the number of folds used for cross-fitting. Default is 5.
#' @param k.cv number of folds used for cross-validation to determine . Default is 10.
#' @param lambda_s fixed step size \eqn{\lambda_s} that is sued for \eqn{s}-function steps in sandwich boosting when \code{variable_steps=FALSE}. Default is 1
#' @param lambda_theta step size for \eqn{\theta}-gradient descent. Default is 0.1.
#' @param m_stop maximum number of sandwich boosting iterations. The number of boosting iterations would therefore in practice be the minimum of value calculated by CV and \code{m_stop}. Default is 100
#' @param init_CCF initalisation of weight class used for \eqn{(s,\theta)}-boosting. Default is \eqn{s=1} and \eqn{\theta=0} (\code{init_CCF = NULL}). Alternatives include "\code{mem_init}" to initialise at an (intercept only) mixed effects model.
#' @return Plots the CV curve for the inputted \code{(lambda_s,m_stop)} pair, alongside a list containing:
#'   \describe{
#'     \item{\code{cv.critical.vals}}{A vector containing the CV criterion (sandwich loss) evaluated at each \code{m_stop} value.}
#'     \item{\code{cv.m_stop}}{The m_stop that minimises the CV criterion (sandwich loss) over the inputted range.}
#'   }
#' @importFrom lme4 lmer
#' @importFrom caret groupKFold
#' @importFrom expm sqrtm
#' @importFrom graphics abline
#' @export
cross.validation.plot <- function(l_formula, l_learner, m_formula, m_learner, s_formula, s_learner, proxyCCF="equicorr", data, K=5L, k.cv=10L, lambda_s=1, lambda_theta=0.1, m_stop=200, init_CCF=NULL) {

  # Reorder rows so that groups are in order
  if (!is.data.frame(data)) {
    stop("Error: Input data must be a dataframe.")
  }
  if (! 'id' %in% colnames(data)) {
    stop("Error: Dataframe data must contain a column labelled 'id' that encodes groups")
  }
  if (proxyCCF=="nested") {
    if (! 'subid' %in% colnames(data)) {
      stop("Error: Dataframe data must contain: a column labelled 'id' that encodes outer groups; a column labelled 'subid' that encodes inner subgroups")
    }
  }
  data <- data[order(data$id),]
  data.groups <- unique(data[,"id"])
  if (proxyCCF == "nested") {
    for (ids in data.groups) {
      sub.data.train <- data[data$id == ids,]
      data[data$id == ids,] <- sub.data.train[order(sub.data.train$subid),]
    }
  }

  cv_folds <- groupKFold(data$id,k=K)#Split into K folds (then can paralellise)
  beta_hat_num_k <- numeric(K); beta_hat_den_k <- numeric(K); V_hat_num_k <- numeric(K)
  for (k in c(1)) {
    cv_fold <- cv_folds[[k]]
    data.nuisance <- data[cv_fold,]
    data.beta <- data[-cv_fold,]

    fit.l <- regress.l(l_formula, l_learner, data.nuisance, data.beta)
    l_residuals.nuisance <- fit.l$l_residuals.nuisance; l_residuals.beta <- fit.l$l_residuals.beta

    fit.m <- regress.m(m_formula, m_learner, data.nuisance, data.beta)
    m_residuals.nuisance <- fit.m$m_residuals.nuisance; m_residuals.beta <- fit.m$m_residuals.beta

    xi_hat.nuisance <- m_residuals.nuisance
    beta_unweighted.nuisance <- sum(l_residuals.nuisance*xi_hat.nuisance)/sum(xi_hat.nuisance^2)
    epsilon_hat.nuisance <- l_residuals.nuisance - beta_unweighted.nuisance*xi_hat.nuisance

    xi_hat.beta <- m_residuals.beta
    beta_unweighted.beta <- sum(l_residuals.beta*xi_hat.beta)/sum(xi_hat.beta^2)
    epsilon_hat.beta <- l_residuals.beta - beta_unweighted.beta*xi_hat.beta

    # Prepare (s,theta) boosting
    boostdf.nuisance <- cbind(data.nuisance, Y_minus_l_hat=l_residuals.nuisance, epsilon_hat=epsilon_hat.nuisance , xi_hat=xi_hat.nuisance)
    boostdf.beta <- cbind(data.beta, Y_minus_l_hat=l_residuals.beta, epsilon_hat=epsilon_hat.beta , xi_hat=xi_hat.beta)

    # Initialise boosting
    if (!is.null(init_CCF)) {
      init_dat <- initialise_boost(boostdf.nuisance, boostdf.beta, beta_unweighted.nuisance, m_formula, init_CCF)
      boostdf.nuisance <- init_dat$boostdf.nuisance
      boostdf.beta <- init_dat$boostdf.beta
    }

    s.boost.output <- s.boost.cv.plots(boostdf.nuisance=boostdf.nuisance, boostdf.beta=boostdf.beta, s_formula=s_formula, s_learner=s_learner, proxyCCF=proxyCCF, m_stop=m_stop, lambda_s=lambda_s, lambda_theta=lambda_theta, k.cv=k.cv)
    cv.crit <- s.boost.output$cv.crit; cv.m_stop <- s.boost.output$cv.m_stop
    plot(seq_len(m_stop), cv.crit, xlab="m_stop", ylab="Scaled CV Criterion (Sandwich Loss)")
    abline(v=cv.m_stop, lty=2)

    if ( (cv.m_stop == m_stop) & (abs(cv.crit[cv.m_stop] - cv.crit[m_stop-1]) > 1e-03) ) warning("m_stop that minimses CV criterion (sandwich loss) is the maximal chosen value. Consider either increasing m_stop or choosing a larger value of lambda_s.")
    if (cv.m_stop <= 10) warning("Sandwich boosting converges very quickly. Consider choosing a smaller value of lambda_s (and / or lambda_theta).")

  }
  return(list(cv.critical.vals=cv.crit, cv.m_stop=cv.m_stop))
}
