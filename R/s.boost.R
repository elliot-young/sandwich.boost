# Functions to carry out (s,\theta)-boosting
#' @importFrom caret groupKFold

s.boost.varstep <- function(boostdf.nuisance, boostdf.beta, s_formula, s_learner, proxyCCF, m_stop, lambda_theta, Lambda_s, mu_s) {
  cv.train.loss <- numeric(m_stop); cv.test.loss <- numeric(m_stop)
  cv_f <- groupKFold(boostdf.nuisance$id,k=10)$Fold1#90% train, 10% test
  boostdf.nuisance.train <- boostdf.nuisance[cv_f,]
  boostdf.nuisance.test <- boostdf.nuisance[-cv_f,]

  #Grouping Details

  if (proxyCCF == "equicorr") {
    grouping.metadata <- grouping.metadata_nonest
    ngradient_cpp <- ngradient_equicorr_cpp
    risk_cpp <- risk_equicorr_cpp
    optimal_step_size_cpp <- optimal_step_size_equicorr_cpp
    theta <- 0
    theta_upper_bound <- 19
  } else if (proxyCCF == "autoreg") {
    grouping.metadata <- grouping.metadata_nonest
    ngradient_cpp <- ngradient_AR_cpp
    risk_cpp <- risk_AR_cpp
    theta <- 0
    theta_upper_bound <- 0.95
  } else if (proxyCCF == "nested") {
    grouping.metadata <- grouping.metadata_nested
    ngradient_cpp <- ngradient_nested_cpp
    risk_cpp <- risk_nested_cpp
    theta <- rep(0,2)
    theta_upper_bound <- 19
  } else {
    stop("Error: proxyCCF should be specified as a valid proxyCCF - equicorr, autoreg or nested")
  }

  grouping.nuisance.train <- grouping.metadata(boostdf.nuisance.train)
  grouping.nuisance.test <- grouping.metadata(boostdf.nuisance.test)
  grouping.beta <- grouping.metadata(boostdf.beta)

  s.train <- rep(1,grouping.nuisance.train$n_obs.metadata)
  s.test <- rep(1,grouping.nuisance.test$n_obs.metadata)
  s.beta <- rep(1,grouping.beta$n_obs.metadata)

  U_S <- all.vars(s_formula)[1]
  if (U_S != "u_s") {
    stop("Error: s_formula should have response in the form 'u_s ~ predictors'")
  }
  X <- all.vars(s_formula)[-1]


  for (m in seq_len(m_stop)) {

    scores <- ngradient_cpp(boostdf.nuisance.train[,"epsilon_hat"], boostdf.nuisance.train[,"xi_hat"], grouping.nuisance.train$I.metadata, grouping.nuisance.train$n_obs.metadata, grouping.nuisance.train$n_i.metadata, s.train, theta)
    u_s <- scores$s
    fit.s <- regress.s(s_formula, s_learner, boostdf.nuisance.train, boostdf.nuisance.test, boostdf.beta, u_s, U_S, X)
    u_s_hat.train <- fit.s$u_s_hat.train; u_s_hat.test <- fit.s$u_s_hat.test; u_s_hat.beta <- fit.s$u_s_hat.beta

    # Line Search + s-boosting
    if (proxyCCF=="equicorr") {
      lambda_optimal <- optimal_step_size_cpp(boostdf.nuisance.train[,"epsilon_hat"], boostdf.nuisance.train[,"xi_hat"], grouping.nuisance.train$I.metadata, grouping.nuisance.train$n_obs.metadata, grouping.nuisance.train$n_i.metadata, s.train, u_s_hat.train, theta)
      lambda_optimal <- max(min(lambda_optimal,Lambda_s[2]), Lambda_s[1])
    } else {
      stop("Error: variable step size only currently supported for proxyCCF = equicorr")
    }

    s.test <- pmax( s.test - mu_s * lambda_optimal * u_s_hat.test, 0.01)
    s.beta <- pmax( s.beta - mu_s * lambda_optimal * u_s_hat.beta, 0.01)
    s.train <- pmax( s.train - mu_s * lambda_optimal * u_s_hat.train, 0.01)
    # Scale to mean 1
    scaling <- max(0.01, mean(c(s.test,s.beta,s.train)))
    s.test <- s.test/scaling
    s.beta <- s.beta/scaling
    s.train <- s.train/scaling

    #theta gradient descent
    u_theta <- scores$theta
    theta <- theta - lambda_theta * u_theta#edit
    theta <- pmin(pmax(theta, 0),theta_upper_bound)

    cv.train.loss[m] <- risk_cpp(boostdf.nuisance.train[,"epsilon_hat"], boostdf.nuisance.train[,"xi_hat"], grouping.nuisance.train$I.metadata, grouping.nuisance.train$n_obs.metadata, grouping.nuisance.train$n_i.metadata, s.train, theta)
    cv.test.loss[m] <- risk_cpp(boostdf.nuisance.test[,"epsilon_hat"], boostdf.nuisance.test[,"xi_hat"], grouping.nuisance.test$I.metadata, grouping.nuisance.test$n_obs.metadata, grouping.nuisance.test$n_i.metadata, s.test, theta)

#    if (m>1) {
#      if ((cv.test.loss[m-1]-cv.test.loss[m])/cv.test.loss[m-1]*100<0) {break}
#    }
    if (m>50) {
      if (1-mean(cv.test.loss[(m-20):m])/mean(cv.test.loss[(m-40):(m-20)])<0 ) {break}
    }

  }

  return(list(s.beta=s.beta, theta=theta))
}

s.boost.constep <- function(boostdf.nuisance, boostdf.beta, s_formula, s_learner, proxyCCF, m_stop, lambda_s, lambda_theta) {
  boostdf.nuisance.train=boostdf.nuisance # to edit
  boostdf.nuisance.test=boostdf.nuisance[1:3,] # to edit

  #Grouping Details

  if (proxyCCF == "equicorr") {
    grouping.metadata <- grouping.metadata_nonest
    ngradient_cpp <- ngradient_equicorr_cpp
    risk_cpp <- risk_equicorr_cpp
    theta <- 0
    theta_upper_bound <- 19
  } else if (proxyCCF == "autoreg") {
    grouping.metadata <- grouping.metadata_nonest
    ngradient_cpp <- ngradient_AR_cpp
    risk_cpp <- risk_AR_cpp
    theta <- 0
    theta_upper_bound <- 0.95
  } else if (proxyCCF == "nested") {
    grouping.metadata <- grouping.metadata_nested
    ngradient_cpp <- ngradient_nested_cpp
    risk_cpp <- risk_nested_cpp
    theta <- rep(0,2)
    theta_upper_bound <- 19
  } else {
    stop("Error: proxyCCF should be specified as a valid proxyCCF - equicorr, autoreg or nested")
  }

  grouping.nuisance.train <- grouping.metadata(boostdf.nuisance.train)
  grouping.nuisance.test <- grouping.metadata(boostdf.nuisance.test)
  grouping.beta <- grouping.metadata(boostdf.beta)

  s.train <- rep(1,grouping.nuisance.train$n_obs.metadata)
  s.test <- rep(1,grouping.nuisance.test$n_obs.metadata)
  s.beta <- rep(1,grouping.beta$n_obs.metadata)

  U_S <- all.vars(s_formula)[1]
#  if (U_S != "u_s") {
#    stop("Error: s_formula should have response in the form 'u_s ~ predictors'")
#  }
  X <- all.vars(s_formula)[-1]


  for (m in seq_len(m_stop)) {

    scores <- ngradient_cpp(boostdf.nuisance.train[,"epsilon_hat"], boostdf.nuisance.train[,"xi_hat"], grouping.nuisance.train$I.metadata, grouping.nuisance.train$n_obs.metadata, grouping.nuisance.train$n_i.metadata, s.train, theta)
    u_s <- scores$s
    fit.s <- regress.s(s_formula, s_learner, boostdf.nuisance.train, boostdf.nuisance.test, boostdf.beta, u_s, U_S, X)
    u_s_hat.train <- fit.s$u_s_hat.train; u_s_hat.test <- fit.s$u_s_hat.test; u_s_hat.beta <- fit.s$u_s_hat.beta

    # s-boosting
    s.test <- pmax( s.test - lambda_s * u_s_hat.test, 0.01)
    s.beta <- pmax( s.beta - lambda_s * u_s_hat.beta, 0.01)
    s.train <- pmax( s.train - lambda_s * u_s_hat.train, 0.01)
    # Scale to mean 1
    scaling <- max(0.01, mean(c(s.test,s.beta,s.train)))
    s.test <- s.test/scaling
    s.beta <- s.beta/scaling
    s.train <- s.train/scaling

    #theta gradient descent
    u_theta <- scores$theta
    theta <- theta - lambda_theta * u_theta#edit
    theta <- pmin(pmax(theta, 0),theta_upper_bound)

  }

  return(list(s.beta=s.beta, theta=theta))
}

s.boost.cv.plots <- function(boostdf.nuisance, boostdf.beta, s_formula, s_learner, m_stop=100, proxyCCF=proxyCCF, early_stopping=TRUE, lambda_s=1, lambda_theta=0.1, k.cv=10) {

  cv.crit <- numeric(m_stop)
  cv_folds <- groupKFold(boostdf.nuisance$id,k=k.cv)

  for (K.cv in 1:k.cv) {
    cv_f <- cv_folds[[K.cv]]
    boostdf.nuisance.train <- boostdf.nuisance[cv_f,]
    boostdf.nuisance.test <- boostdf.nuisance[-cv_f,]

    # Grouping metadata
    if (proxyCCF == "equicorr") {
      grouping.metadata <- grouping.metadata_nonest
      ngradient_cpp <- ngradient_equicorr_cpp
      risk_cpp <- risk_equicorr_cpp
      theta <- 0
      theta_upper_bound <- 19
    } else if (proxyCCF == "autoreg") {
      grouping.metadata <- grouping.metadata_nonest
      ngradient_cpp <- ngradient_AR_cpp
      risk_cpp <- risk_AR_cpp
      theta <- 0
      theta_upper_bound <- 0.95
    } else if (proxyCCF == "nested") {
      grouping.metadata <- grouping.metadata_nested
      ngradient_cpp <- ngradient_nested_cpp
      risk_cpp <- risk_nested_cpp
      theta <- rep(0,2)
      theta_upper_bound <- 19
    } else {
      stop("Error: proxyCCF should be specified as a valid proxyCCF - equicorr, autoreg or nested")
    }

    #Grouping Details
    grouping.nuisance.train <- grouping.metadata(boostdf.nuisance.train)
    groups.nuisance.train <- grouping.nuisance.train$groups.metadata
    I.nuisance.train <- grouping.nuisance.train$I.metadata
    n_i.nuisance.train <- grouping.nuisance.train$n_i.metadata
    n_obs.nuisance.train <- grouping.nuisance.train$n_obs.metadata

    grouping.nuisance.test <- grouping.metadata(boostdf.nuisance.test)
    groups.nuisance.test <- grouping.nuisance.test$groups.metadata
    I.nuisance.test <- grouping.nuisance.test$I.metadata
    n_i.nuisance.test <- grouping.nuisance.test$n_i.metadata
    n_obs.nuisance.test <- grouping.nuisance.test$n_obs.metadata

    grouping.beta <- grouping.metadata(boostdf.beta)
    groups.beta <- grouping.beta$groups.metadata
    I.beta <- grouping.beta$I.metadata
    n_i.beta <- grouping.beta$n_i.metadata
    n_obs.beta <- grouping.beta$n_obs.metadata

    #Boost

    s.train <- rep(1,n_obs.nuisance.train)
    #theta <- 1/mean(n_i.nuisance.train) #CHANGE TO ALLOW USER INPUT!!!
    theta <- 0
    s.test <- rep(1,n_obs.nuisance.test)
    s.beta <- rep(1,n_obs.beta)

    U_S <- all.vars(s_formula)[1]
    if (U_S != "u_s") {
      stop("Error: s_formula should have response in the form 'u_s ~ predictors'")
    }
    X <- all.vars(s_formula)[-1]

    for (m in seq_len(m_stop)) {
      scores <- ngradient_cpp(boostdf.nuisance.train[,"epsilon_hat"], boostdf.nuisance.train[,"xi_hat"], I.nuisance.train, n_obs.nuisance.train, n_i.nuisance.train, s.train, theta)
      u_s <- scores$s
      fit.s <- regress.s(s_formula, s_learner, boostdf.nuisance.train, boostdf.nuisance.test, boostdf.beta, u_s, U_S, X)
      u_s_hat.train <- fit.s$u_s_hat.train; u_s_hat.test <- fit.s$u_s_hat.test; u_s_hat.beta <- fit.s$u_s_hat.beta

      # s-boosting
      s.test <- pmax( s.test - lambda_s * u_s_hat.test, 0.01)
      s.beta <- pmax( s.beta - lambda_s * u_s_hat.beta, 0.01)
      s.train <- pmax( s.train - lambda_s * u_s_hat.train, 0.01)
      # Scale to mean 1
      scaling <- max(0.01, mean(c(s.test,s.beta,s.train)))
      s.test <- s.test/scaling
      s.beta <- s.beta/scaling
      s.train <- s.train/scaling

      #theta gradient descent
      u_theta <- scores$theta
      theta <- theta - lambda_theta * u_theta#edit
      theta <- pmin(pmax(theta, 0),theta_upper_bound)

      cv.crit[m] <- cv.crit[m] + risk_cpp(boostdf.nuisance.test[,"epsilon_hat"], boostdf.nuisance.test[,"xi_hat"], I.nuisance.test, n_obs.nuisance.test, n_i.nuisance.test, s.test, theta)

    }

  }

  #plot(seq_len(m_stop), cv.crit, xlab="m_stop", ylab="CV Criterion")

  cv.m_stop <- which.min(cv.crit)

  return(list(cv.crit=cv.crit, cv.m_stop=cv.m_stop))
}
