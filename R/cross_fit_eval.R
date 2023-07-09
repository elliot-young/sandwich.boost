# Functions to evaluate risk through R (as opposed to C++ code).
# Redundant - For historical/checking purposes only.

cross_fit_eval <- function(boostdf.beta, s.beta, theta, proxyCCF) {
  if (proxyCCF == "equicorr") {
    groups.beta <- unique(boostdf.beta[,"id"])
    I.beta <- length(groups.beta)
    beta_num = 0; beta_den = 0; V_num = 0
    for (i in seq_len(I.beta)) {
      Y_minus_l_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"Y_minus_l_hat"]
      xi_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"xi_hat"]
      epsilon_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"epsilon_hat"]
      s_i <- s.beta[is.element(boostdf.beta$id,groups.beta[i])]
      N <- length(s_i)
      A1 <- sum(s_i^2*xi_hat_i^2) - (theta/(1+N*theta))*(sum(s_i*xi_hat_i))^2
      A2_Y_minus_l_hat <- sum(s_i^2*Y_minus_l_hat_i*xi_hat_i) - (theta/(1+N*theta))*sum(s_i*Y_minus_l_hat_i)*sum(s_i*xi_hat_i)
      A2 <- sum(s_i^2*epsilon_hat_i*xi_hat_i) - (theta/(1+N*theta))*sum(s_i*epsilon_hat_i)*sum(s_i*xi_hat_i)
      beta_num = beta_num + A2_Y_minus_l_hat
      beta_den = beta_den + A1
      V_num = V_num + A2^2
    }
  } else if (proxyCCF == "autoreg") {
    groups.beta <- unique(boostdf.beta[,"id"])
    I.beta <- length(groups.beta)
    beta_num = 0; beta_den = 0; V_num = 0
    for (i in seq_len(I.beta)) {
      Y_minus_l_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"Y_minus_l_hat"]
      xi_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"xi_hat"]
      epsilon_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"epsilon_hat"]
      s_i <- s.beta[is.element(boostdf.beta$id,groups.beta[i])]
      N <- length(s_i)
      if (N>2) {
        beta1 <- sum((s_i[2:(N-1)]*xi_hat_i[2:(N-1)])^2)
        alpha1 <- beta1 + (s_i[1]*xi_hat_i[1])^2 + (s_i[N]*xi_hat_i[N])^2
        gamma1 <- 2*sum(s_i[2:N]*s_i[1:(N-1)]*xi_hat_i[2:N]*xi_hat_i[1:(N-1)])
        beta2_Y_minus_l_hat <- sum(s_i[2:(N-1)]^2*xi_hat_i[2:(N-1)]*Y_minus_l_hat_i[2:(N-1)])
        alpha2_Y_minus_l_hat <- beta2_Y_minus_l_hat + s_i[1]^2*xi_hat_i[1]*Y_minus_l_hat_i[1] + s_i[N]^2*xi_hat_i[N]*Y_minus_l_hat_i[N]
        gamma2_Y_minus_l_hat <- sum(s_i[2:N]*s_i[1:(N-1)]*(xi_hat_i[2:N]*Y_minus_l_hat_i[1:(N-1)] + Y_minus_l_hat_i[2:N]*xi_hat_i[1:(N-1)]))
        beta2 <- sum(s_i[2:(N-1)]^2*xi_hat_i[2:(N-1)]*epsilon_hat_i[2:(N-1)])
        alpha2 <- beta2 + s_i[1]^2*xi_hat_i[1]*epsilon_hat_i[1] + s_i[N]^2*xi_hat_i[N]*epsilon_hat_i[N]
        gamma2 <- sum(s_i[2:N]*s_i[1:(N-1)]*(xi_hat_i[2:N]*epsilon_hat_i[1:(N-1)] + epsilon_hat_i[2:N]*xi_hat_i[1:(N-1)]))
        A1 <- alpha1 + theta^2*beta1 - theta*gamma1
        A2_Y_minus_l_hat <- alpha2_Y_minus_l_hat + theta^2*beta2_Y_minus_l_hat - theta*gamma2_Y_minus_l_hat
        A2 <- alpha2 + theta^2*beta2 - theta*gamma2
      } else if (N==2) {
        alpha1 <- (s_i[1]*xi_hat_i[1])^2 + (s_i[N]*xi_hat_i[N])^2
        gamma1 <- 2*sum(s_i[N]*s_i[N-1]*xi_hat_i[N]*xi_hat_i[N-1])
        alpha2_Y_minus_l_hat <- s_i[1]^2*xi_hat_i[1]*Y_minus_l_hat_i[1] + s_i[N]^2*xi_hat_i[N]*Y_minus_l_hat_i[N]
        gamma2_Y_minus_l_hat <- sum(s_i[N]*s_i[N-1]*(xi_hat_i[N]*Y_minus_l_hat_i[N-1] + Y_minus_l_hat_i[N]*xi_hat_i[N-1]))
        alpha2 <- s_i[1]^2*xi_hat_i[1]*epsilon_hat_i[1] + s_i[N]^2*xi_hat_i[N]*epsilon_hat_i[N]
        gamma2 <- sum(s_i[N]*s_i[N-1]*(xi_hat_i[N]*epsilon_hat_i[N-1] + epsilon_hat_i[N]*xi_hat_i[N-1]))
        A1 <- alpha1 - theta*gamma1
        A2_Y_minus_l_hat <- alpha2_Y_minus_l_hat - theta*gamma2_Y_minus_l_hat
        A2 <- alpha2 - theta*gamma2
      } else if (N==1) {
        alpha1 <- (s_i[1]*xi_hat_i[1])^2 + (s_i[N]*xi_hat_i[N])^2
        alpha2_Y_minus_l_hat <- s_i[1]^2*xi_hat_i[1]*Y_minus_l_hat_i[1] + s_i[N]^2*xi_hat_i[N]*Y_minus_l_hat_i[N]
        alpha2 <- s_i[1]^2*xi_hat_i[1]*epsilon_hat_i[1] + s_i[N]^2*xi_hat_i[N]*epsilon_hat_i[N]
        A1 <- alpha1
        A2_Y_minus_l_hat <- alpha2_Y_minus_l_hat
        A2 <- alpha2
      }
      beta_num = beta_num + A2_Y_minus_l_hat
      beta_den = beta_den + A1
      V_num = V_num + A2^2
    }
  } else if (proxyCCF == "nested"){
    groups.beta <- unique(boostdf.beta[,"id"])
    I.beta <- length(groups.beta)
    beta_num = 0; beta_den = 0; V_num = 0
    for (i in seq_len(I.beta)) {
      sub.boostdf.beta <- boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),]
      sub.s.beta <- s.beta[is.element(boostdf.beta$id,groups.beta[i])]
      sub.groups.beta <- unique(sub.boostdf.beta[,"subid"])
      J.beta <- length(sub.groups.beta)
      b1 <- b2 <- b3 <- c1 <- c2 <- c3 <- n_ij <- rep(0,J.beta)
      for (j in seq_len(J.beta)) {
        Y_minus_l_hat_ij = sub.boostdf.beta[is.element(sub.boostdf.beta$subid,sub.groups.beta[j]),"Y_minus_l_hat"]
        xi_hat_ij = sub.boostdf.beta[is.element(sub.boostdf.beta$subid,sub.groups.beta[j]),"xi_hat"]
        epsilon_hat_ij = sub.boostdf.beta[is.element(sub.boostdf.beta$subid,sub.groups.beta[j]),"epsilon_hat"]
        s_ij = sub.s.beta[is.element(sub.boostdf.beta$subid,sub.groups.beta[j])]
        b1[j] = sum(s_ij*xi_hat_ij)
        b2[j] = sum(s_ij*epsilon_hat_ij)
        b3[j] = sum(s_ij*Y_minus_l_hat_ij)
        n_ij[j] = length(s_ij)
        c1[j] = sum(s_ij^2*xi_hat_ij^2)
        c2[j] = sum(s_ij^2*xi_hat_ij*epsilon_hat_ij)
        c3[j] = sum(s_ij^2*xi_hat_ij*Y_minus_l_hat_ij)
      }
      phi <- sum(n_ij/(1+theta[1]*n_ij))
      A1 <- (1+theta[1]+theta[2])*sum(c1) - sum(theta[1]/(1+theta[1]*n_ij)*b1^2) - theta[2]/(1+theta[2]*phi)*(sum(b1/(1+theta[1]*n_ij)))^2
      A2 <- (1+theta[1]+theta[2])*sum(c2) - sum(theta[1]/(1+theta[1]*n_ij)*b1*b2) - theta[2]/(1+theta[2]*phi)*sum(b1/(1+theta[1]*n_ij))*sum(b2/(1+theta[1]*n_ij))
      A2_Y_minus_l_hat <- (1+theta[1]+theta[2])*sum(c3) - sum(theta[1]/(1+theta[1]*n_ij)*b1*b3) - theta[2]/(1+theta[2]*phi)*sum(b1/(1+theta[1]*n_ij))*sum(b3/(1+theta[1]*n_ij))
      beta_num = beta_num + A2_Y_minus_l_hat
      beta_den = beta_den + A1
      V_num = V_num + A2^2
    }
  }
  return(list(beta_num=beta_num, beta_den=beta_den, V_num=V_num))
}
