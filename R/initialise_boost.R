# Function to initialise sandwich boosting by transforming residuals
#' @importFrom lme4 lmer
#' @importFrom expm sqrtm

initialise_boost <- function(df.nuisance, df.beta, beta.unw.nuis, m_formula, init_CCF) {

  if (init_CCF == "mem") {
    lmerD <- all.vars(m_formula)[1];   lmerX <- all.vars(m_formula)[-1]
    lmer_response <- df.nuisance$Y_minus_l_hat - beta.unw.nuis*(df.nuisance$xi_hat-df.nuisance[,lmerD])
    data.lmer <- data.frame(lmer_D=df.nuisance[,lmerD], lmer_response=lmer_response, id=df.nuisance[,"id"])
    lmer_fit <- lmer(lmer_response ~ lmer_D - 1 + (1|id), data=data.lmer)
    lmer_sigma <- sigma(lmer_fit)
    lmer_cov_matrix <- summary(lmer_fit)$varcor$id
    groups.nuisance <- unique(df.nuisance[,"id"])
    for (i in groups.nuisance) {
      x = df.nuisance[df.nuisance$id==i,lmerX]
      lengthX <- length(lmerX); N <- length(x)/lengthX
      Sigma <- matrix(c(rep(1,N)),nrow=N,ncol=1) %*% lmer_cov_matrix %*% t(matrix(c(rep(1,N)),nrow=N,ncol=1)) + lmer_sigma*diag(N)
      MEM_trans <- sqrtm(solve(Sigma))
      df.nuisance[is.element(df.nuisance$id,i),"Y_minus_l_hat"] = MEM_trans %*% df.nuisance[is.element(df.nuisance$id,i),"Y_minus_l_hat"]
      df.nuisance[is.element(df.nuisance$id,i),"epsilon_hat"] = MEM_trans %*% df.nuisance[is.element(df.nuisance$id,i),"epsilon_hat"]
      df.nuisance[is.element(df.nuisance$id,i),"xi_hat"] = MEM_trans %*% df.nuisance[is.element(df.nuisance$id,i),"xi_hat"]
    }
    groups.beta <- unique(df.beta[,"id"])
    for (i in groups.beta) {
      x = df.beta[df.beta$id==i,lmerX]
      lengthX <- length(lmerX); N <- length(x)/lengthX
      Sigma <- matrix(c(rep(1,N)),nrow=N,ncol=1) %*% lmer_cov_matrix %*% t(matrix(c(rep(1,N)),nrow=N,ncol=1)) + lmer_sigma*diag(N)
      MEM_trans <- sqrtm(solve(Sigma))
      df.beta[is.element(df.beta$id,i),"Y_minus_l_hat"] = MEM_trans %*% df.beta[is.element(df.beta$id,i),"Y_minus_l_hat"]
      df.beta[is.element(df.beta$id,i),"epsilon_hat"] = MEM_trans %*% df.beta[is.element(df.beta$id,i),"epsilon_hat"]
      df.beta[is.element(df.beta$id,i),"xi_hat"] = MEM_trans %*% df.beta[is.element(df.beta$id,i),"xi_hat"]
    }
  } else if (init_CCF == "mem_int") {
    lmerD <- all.vars(m_formula)[1];   lmerX <- all.vars(m_formula)[-1]
    lmer_response <- df.nuisance$Y_minus_l_hat - beta.unw.nuis*(df.nuisance$xi_hat - df.nuisance[,lmerD])
    data.lmer <- data.frame(lmer_D=df.nuisance[,lmerD], lmer_response=lmer_response, id=df.nuisance[,"id"])
    lmer_fit <- lmer(lmer_response ~ lmer_D - 1 + (1|id), data=data.lmer)
    lmer_sigma <- sigma(lmer_fit)
    lmer_cov_matrix <- summary(lmer_fit)$varcor$id
    groups.nuisance <- unique(df.nuisance[,"id"])
    for (i in groups.nuisance) {
      x = df.nuisance[df.nuisance$id==i,lmerX]
      lengthX <- length(lmerX); N <- length(x)/lengthX
      Sigma <- matrix(c(rep(1,N)),nrow=N,ncol=1) %*% lmer_cov_matrix %*% t(matrix(c(rep(1,N)),nrow=N,ncol=1)) + lmer_sigma*diag(N)
      MEM_trans <- sqrtm(solve(Sigma))
      df.nuisance[is.element(df.nuisance$id,i),"Y_minus_l_hat"] = MEM_trans %*% df.nuisance[is.element(df.nuisance$id,i),"Y_minus_l_hat"]
      df.nuisance[is.element(df.nuisance$id,i),"epsilon_hat"] = MEM_trans %*% df.nuisance[is.element(df.nuisance$id,i),"epsilon_hat"]
      df.nuisance[is.element(df.nuisance$id,i),"xi_hat"] = MEM_trans %*% df.nuisance[is.element(df.nuisance$id,i),"xi_hat"]
    }
    groups.beta <- unique(df.beta[,"id"])
    for (i in groups.beta) {
      x = df.beta[df.beta$id==i,lmerX]
      lengthX <- length(lmerX); N <- length(x)/lengthX
      Sigma <- matrix(c(rep(1,N)),nrow=N,ncol=1) %*% lmer_cov_matrix %*% t(matrix(c(rep(1,N)),nrow=N,ncol=1)) + lmer_sigma*diag(N)
      MEM_trans <- sqrtm(solve(Sigma))
      df.beta[is.element(df.beta$id,i),"Y_minus_l_hat"] = MEM_trans %*% df.beta[is.element(df.beta$id,i),"Y_minus_l_hat"]
      df.beta[is.element(df.beta$id,i),"epsilon_hat"] = MEM_trans %*% df.beta[is.element(df.beta$id,i),"epsilon_hat"]
      df.beta[is.element(df.beta$id,i),"xi_hat"] = MEM_trans %*% df.beta[is.element(df.beta$id,i),"xi_hat"]
    }
  }

  return(return(list(boostdf.nuisance=df.nuisance, boostdf.beta=df.beta)))
}
