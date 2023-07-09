#' Sandwich Boosting Estimator for the PLR Model
#'
#' @param l_formula a two-sided formula object describing the regression model for \eqn{E[Y|X=x]}.
#' @param l_learner a string specifying the regression method to fit the regression of \eqn{Y} on \eqn{X} as given by \code{l_formula}. Notable options are \code{gam} and \code{randomforest}
#' @param m_formula a two-sided formula object describing the regression model for \eqn{E[D|X=x]}.
#' @param m_learner a string specifying the regression method to fit the regression of \eqn{D} on \eqn{X} as given by \code{m_formula}.
#' @param s_formula a two-sided formula object describing the regression model for a single boosting iteration for \eqn{s}-function sandwich boosting.
#' @param s_learner a string specifying the regression method to fit the a single boosting iteration for \eqn{s}-function sandwich boosting as given by \code{l_formula}.
#' @param proxyCCF a string spescifying which working correlation parametrisation to use. Must take the form "\code{equicorr}", "\code{autoreg}" or "\code{nested}".
#' @param data a data frame containing the variables for the grouped PLR model. Note that group ID must be given by the column \code{id}.
#' @param K the number of folds used for cross-fitting. Default is 5.
#' @param S the number of repeats to mitigate the randomness in the estimator on the sample splits used for cross-fitting. Default is 5.
#' @param variable_steps whether sandwich boosting on the \eqn{s}-function is performed with fixed step size \eqn{\lambda_s} or variable step size. Default is \code{FALSE}.
#' @param m_stop maximum number of sandwich boosting iterations. The number of boosting iterations would therefore in practice be the minimum of value calculated by CV and \code{m_stop}. Default is 100
#' @param lambda_s fixed step size \eqn{\lambda_s} that is sued for \eqn{s}-function steps in sandwich boosting when \code{variable_steps=FALSE}. Default is 1
#' @param lambda_theta step size for \eqn{\theta}-gradient descent. Default is 0.1.
#' @param Lambda_s a closed set of positive step sizes to restrict step sizes for \eqn{s}-descent when \code{variable_steps=TRUE}. Default is \eqn{[0.01,10]}
#' @param mu_s shrinkage parameter \eqn{\mu_s\in[0,1]} for \eqn{s}-function sandwich boosting step. Default is 0.1.
#' @param init_CCF initalisation of weight class used for \eqn{(s,\theta)}-boosting. Default is \eqn{s=1} and \eqn{\theta=0} (\code{init_CCF = NULL}). Alternatives include "\code{mem}" and "\code{mem_int}" for mixed effects model .
#'
#' @return A list containing:
#'   \describe{
#'     \item{\code{beta_hat}}{The estimator for \eqn{\beta}.}
#'     \item{\code{std_err}}{Huber robust estimate of the standard error of the \eqn{\beta}-estimator.}
#'   }
#' @importFrom lme4 lmer
#' @importFrom caret groupKFold
#' @importFrom expm sqrtm
#' @importFrom stats median
#' @export
PLR_sandboost <- function(l_formula, l_learner, m_formula, m_learner, s_formula, s_learner, proxyCCF="equicorr", data, K=5L, S=5L, variable_steps = FALSE, m_stop=100, lambda_s=1, lambda_theta=0.1, Lambda_s=c(0.01,10), mu_s=0.1, init_CCF=NULL) {

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

  beta_hat <- V_hat_div_n <- numeric(S)
  for (s in seq_len(S)) {
    cv_folds <- groupKFold(data$id,k=K)#Split into K folds (then can paralellise)
    beta_hat_num_k <- beta_hat_den_k <- V_hat_num_k <- numeric(K)
    for (k in seq_len(K)) {
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

      if (!variable_steps) s.boost.output <- s.boost.constep(boostdf.nuisance=boostdf.nuisance, boostdf.beta=boostdf.beta, s_formula=s_formula, s_learner=s_learner, proxyCCF=proxyCCF, m_stop=m_stop, lambda_s=lambda_s, lambda_theta=lambda_theta)
        else s.boost.output <- s.boost.varstep(boostdf.nuisance=boostdf.nuisance, boostdf.beta=boostdf.beta, s_formula=s_formula, s_learner=s_learner, proxyCCF=proxyCCF, m_stop=m_stop, lambda_theta=lambda_theta, Lambda_s=Lambda_s, mu_s=mu_s)

      #s.boost.output <- s.boost(boostdf.nuisance=boostdf.nuisance, boostdf.beta=boostdf.beta, s_formula=s_formula, s_learner=s_learner, proxyCCF=proxyCCF, variable_steps=variable_steps, m_stop=m_stop, lambda_s=lambda_s, lambda_theta=lambda_theta, Lambda_theta=Lambda_theta, mu_s=mu_s)
      s.beta <- s.boost.output$s.beta; theta <- s.boost.output$theta

      cross_fit_evals <- cross_fit_eval(boostdf.beta=boostdf.beta, s.beta=s.beta, theta=theta, proxyCCF=proxyCCF)
      beta_hat_num_k[k] <- cross_fit_evals$beta_num
      beta_hat_den_k[k] <- cross_fit_evals$beta_den
      V_hat_num_k[k] <- cross_fit_evals$V_num

    }
    beta_hat[s] <- sum(beta_hat_num_k)/sum(beta_hat_den_k)
    V_hat_div_n[s] <- sum(V_hat_num_k)/(sum(beta_hat_den_k))^2
  }

#  PLR_output <- list(beta_hat = median(beta_hat),
#                     std_err = sqrt(median(V_hat_div_n+(beta_hat-median(beta_hat))^2)))

  se <- sqrt(median(V_hat_div_n+(beta_hat-median(beta_hat))^2))
  beta <- median(beta_hat)
  zval <- tval <- beta/se
  tab <- cbind(Estimate = beta,
               StdErr = se,
               z.value = zval,
               p.value = 2*pnorm(-abs(zval)))
  rownames(tab) <- all.vars(m_formula)[1] #c("Treatment (D)")
  colnames(tab) <- c("Estimate", "Std. Err", "z value", "Pr(>|z|)")
  res <- list(coefficients=tab)
  class(res) <- "sandwichboost"
  res

  return(res)
}

