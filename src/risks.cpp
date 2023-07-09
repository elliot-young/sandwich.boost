#include <Rcpp.h>
using namespace Rcpp;


// Calculates risks for sandwich boosting CV

// Risks for equicorrelated proxyCCF
// [[Rcpp::export]]
double risk_equicorr_cpp(NumericVector epsilon_hat, NumericVector xi_hat, int I, int n_obs, NumericVector n_i, NumericVector s, double theta) {
  NumericVector A1(I);
  NumericVector A2(I);
  int n_cumsum_start = 0;
  for(int i = 0; i < I; i++) {
    NumericVector s_i(n_i(i));
    NumericVector xi_hat_i(n_i(i));
    NumericVector epsilon_hat_i(n_i(i));
    for(int gg = 0; gg < n_i(i); gg++) {
      s_i(gg) = s(n_cumsum_start + gg);
      xi_hat_i(gg) = xi_hat(n_cumsum_start + gg);
      epsilon_hat_i(gg) = epsilon_hat(n_cumsum_start + gg);
    }
    A1(i) = sum(pow(s_i,2.0)*pow(xi_hat_i,2.0))-theta/(1+theta*n_i(i))*pow((sum(s_i*xi_hat_i)),2.0);
    A2(i) = sum(pow(s_i,2.0)*xi_hat_i*epsilon_hat_i)-theta/(1+theta*n_i(i))*(sum(s_i*xi_hat_i))*(sum(s_i*epsilon_hat_i));
    n_cumsum_start += n_i(i);
  }

  NumericVector s_scores(n_obs);
  double sumA1 = sum(A1);
  double sumA2sq = sum(pow(A2,2.0));

  double risk = n_obs*sumA2sq/pow(sumA1,2.0);

  return risk;
}

// Risks for autoregressive proxyCCF
// [[Rcpp::export]]
double risk_AR_cpp(NumericVector epsilon_hat, NumericVector xi_hat, int I, int n_obs, NumericVector n_i, NumericVector s, double theta) {
    NumericVector A1(I);
    NumericVector A2(I);
    NumericVector alpha1(I);
    NumericVector beta1(I);
    NumericVector gamma1(I);
    NumericVector alpha2(I);
    NumericVector beta2(I);
    NumericVector gamma2(I);
    int n_cumsum_start = 0;
    for(int i = 0; i < I; i++) {
      NumericVector s_i(n_i(i));
      NumericVector xi_hat_i(n_i(i));
      NumericVector epsilon_hat_i(n_i(i));
      for(int gg = 0; gg < n_i(i); gg++) {
        s_i(gg) = s(n_cumsum_start + gg);
        xi_hat_i(gg) = xi_hat(n_cumsum_start + gg);
        epsilon_hat_i(gg) = epsilon_hat(n_cumsum_start + gg);
      }
      if (n_i(i) > 2){
        NumericVector s_i_sub1(n_i(i)-1);
        NumericVector s_i_sub2(n_i(i)-1);
        NumericVector xi_hat_i_sub1(n_i(i)-1);
        NumericVector xi_hat_i_sub2(n_i(i)-1);
        NumericVector epsilon_hat_i_sub1(n_i(i)-1);
        NumericVector epsilon_hat_i_sub2(n_i(i)-1);
        for(int gg = 0; gg < (n_i(i)-1); gg++) {
          s_i_sub1(gg) = s_i(gg);
          s_i_sub2(gg) = s_i(gg+1);
          xi_hat_i_sub1(gg) = xi_hat_i(gg);
          xi_hat_i_sub2(gg) = xi_hat_i(gg+1);
          epsilon_hat_i_sub1(gg) = epsilon_hat_i(gg);
          epsilon_hat_i_sub2(gg) = epsilon_hat_i(gg+1);
        }
        alpha1(i) = sum(pow(s_i,2.0)*pow(xi_hat_i,2.0));
        beta1(i) = alpha1(i) - pow(s_i(0)*xi_hat_i(0),2.0) - pow(s_i(n_i(i)-1)*xi_hat_i(n_i(i)-1),2.0);
        gamma1(i) = 2*sum(s_i_sub1*s_i_sub2*xi_hat_i_sub1*xi_hat_i_sub2);
        alpha2(i) = sum(pow(s_i,2.0)*xi_hat_i*epsilon_hat_i);
        beta2(i) = alpha2(i) - pow(s_i(0),2.0)*xi_hat_i(0)*epsilon_hat_i(0) - pow(s_i(n_i(i)-1),2.0)*xi_hat_i(n_i(i)-1)*epsilon_hat_i(n_i(i)-1);
        gamma2(i) = sum(s_i_sub1*s_i_sub2*(xi_hat_i_sub2*epsilon_hat_i_sub1+xi_hat_i_sub1*epsilon_hat_i_sub2));
        A1(i) = alpha1(i) + pow(theta,2.0)*beta1(i) - theta*gamma1(i);
        A2(i) = alpha2(i) + pow(theta,2.0)*beta2(i) - theta*gamma2(i);
      }
      else if (n_i(i) == 2){
        alpha1(i) = sum(pow(s_i,2.0)*pow(xi_hat_i,2.0));
        gamma1(i) = 2*(s_i(1)*s_i(0)*xi_hat_i(1)*xi_hat_i(0));
        alpha2(i) = sum(pow(s_i,2.0)*xi_hat_i*epsilon_hat_i);
        gamma2(i) = s_i(1)*s_i(0)*(xi_hat_i(0)*epsilon_hat_i(1)+xi_hat_i(1)*epsilon_hat_i(0));
        A1(i) = alpha1(i) - theta*gamma1(i);
        A2(i) = alpha2(i) - theta*gamma2(i);
      }
      else if (n_i(i) == 1){
        alpha1(i) = 2*pow(s_i(0),2.0)*pow(xi_hat_i(0),2.0);
        alpha2(i) = 2*pow(s_i(0),2.0)*xi_hat_i(0)*epsilon_hat_i(0);
        A1(i) = (1-pow(theta,2.0))*alpha1(i);
        A2(i) = (1-pow(theta,2.0))*alpha2(i);
      }
      n_cumsum_start += n_i(i);
    }

    NumericVector s_scores(n_obs);
    double sumA1 = sum(A1);
    double sumA2sq = sum(pow(A2,2.0));

    double risk = n_obs*sumA2sq/pow(sumA1,2.0);

    return risk;
}

// Risks for nested proxyCCF
// [[Rcpp::export]]
double risk_nested_cpp(NumericVector epsilon_hat, NumericVector xi_hat, NumericVector J, int n_obs, NumericMatrix n_ij, NumericVector s, NumericVector theta) {
  int I = J.length();
  NumericVector phi(I);
  NumericVector phidiff(I);
  NumericVector B1_sum1(I);
  NumericVector B1_sum4(I);
  NumericVector B2_sum1(I);
  NumericVector B2_sum4(I);
  NumericVector A1(I);
  NumericVector A2(I);
  NumericVector C1(I);
  NumericVector C2(I);
  NumericVector A3(n_obs);
  NumericVector A4(n_obs);
  int n_cumsum_start = 0;
  for(int i = 0; i < I; i++) {
    NumericVector b1(J(i));
    NumericVector b2(J(i));
    NumericVector c1(J(i));
    NumericVector c2(J(i));
    NumericVector N_ij(J(i));
    for(int gg = 0; gg < J(i); gg++) {
      N_ij(gg) = n_ij(i,gg);
    }
    for (int j = 0; j < J(i); j++) {
      NumericVector s_i(N_ij(j));
      NumericVector xi_hat_i(N_ij(j));
      NumericVector epsilon_hat_i(N_ij(j));
      for(int gg = 0; gg < N_ij(j); gg++) {
        s_i(gg) = s(n_cumsum_start + gg);
        xi_hat_i(gg) = xi_hat(n_cumsum_start + gg);
        epsilon_hat_i(gg) = epsilon_hat(n_cumsum_start + gg);
      }
      n_cumsum_start += N_ij(j);
      b1(j) = sum(s_i*xi_hat_i);
      b2(j) = sum(s_i*epsilon_hat_i);
      c1(j) = sum(pow(s_i*xi_hat_i,2.0));
      c2(j) = sum(pow(s_i,2.0)*epsilon_hat_i*xi_hat_i);
    }
    phi(i) = sum(N_ij/(1+theta(0)*N_ij));
    phidiff(i) = -sum(pow(N_ij/(1+theta(0)*N_ij),2.0));
    B1_sum1(i) = sum(b1/(1+theta(0)*N_ij));
    B1_sum4(i) = sum(N_ij/(1+theta(0)*N_ij)*pow(b1,2.0));
    B2_sum1(i) = sum(b2/(1+theta(0)*N_ij));
    B2_sum4(i) = sum(N_ij/(1+theta(0)*N_ij)*b1*b2);
    C1(i) = sum(c1);
    C2(i) = sum(c2);
    A1(i) = (1+theta(0)+theta(1))*C1(i) - B1_sum4(i) - theta(1)/(1+theta(1)*phi(i))*pow(B1_sum1(i),2.0);
    A2(i) = (1+theta(0)+theta(1))*C2(i) - B2_sum4(i) - theta(1)/(1+theta(1)*phi(i))*B1_sum1(i)*B2_sum1(i);
  }

  double sumA1 = sum(A1);
  double sumA2sq = sum(pow(A2,2.0));

  double risk = n_obs*sumA2sq/pow(sumA1,2.0);

  return risk;
}


