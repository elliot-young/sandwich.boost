#include <Rcpp.h>
using namespace Rcpp;


// Calculates optimal step size for sandwich boosting

// Optimal step size for equicorrelated proxyCCF
// [[Rcpp::export]]

double optimal_step_size_equicorr_cpp(NumericVector epsilon_hat, NumericVector xi_hat, int I, int n_obs, NumericVector n_i, NumericVector s, NumericVector h, double theta) {
  NumericVector A1(I);
  NumericVector A2(I);
  NumericVector Ab(I);
  NumericVector Ac(I);
  NumericVector Bb(I);
  NumericVector Bc(I);
  int n_cumsum_start = 0;
  for(int i = 0; i < I; i++) {
    NumericVector s_i(n_i(i));
    NumericVector h_i(n_i(i));
    NumericVector xi_hat_i(n_i(i));
    NumericVector epsilon_hat_i(n_i(i));
    for(int gg = 0; gg < n_i(i); gg++) {
      s_i(gg) = s(n_cumsum_start + gg);
      h_i(gg) = h(n_cumsum_start + gg);
      xi_hat_i(gg) = xi_hat(n_cumsum_start + gg);
      epsilon_hat_i(gg) = epsilon_hat(n_cumsum_start + gg);
    }
    A1(i) = sum(pow(s_i,2.0)*pow(xi_hat_i,2.0))-theta/(1+theta*n_i(i))*pow((sum(s_i*xi_hat_i)),2.0);
    A2(i) = sum(pow(s_i,2.0)*xi_hat_i*epsilon_hat_i)-theta/(1+theta*n_i(i))*(sum(s_i*xi_hat_i))*(sum(s_i*epsilon_hat_i));
    Ab(i) = 2*sum(xi_hat_i*epsilon_hat_i*h_i*s_i)-theta/(1+theta*n_i(i))*(sum(s_i*xi_hat_i)*sum(h_i*epsilon_hat_i)+sum(h_i*xi_hat_i)*sum(s_i*epsilon_hat_i));
    Ac(i) = sum(pow(xi_hat_i,2.0)*pow(h_i,2.0))-theta/(1+theta*n_i(i))*pow((sum(h_i*xi_hat_i)),2.0);
    Bb(i) = sum(pow(xi_hat_i,2.0)*h_i*s_i)-theta/(1+theta*n_i(i))*(sum(h_i*xi_hat_i))*(sum(s_i*xi_hat_i));
    Bc(i) = sum(pow(xi_hat_i,2.0)*pow(h_i,2.0))-theta/(1+theta*n_i(i))*pow(sum(h_i*xi_hat_i),2.0);
    n_cumsum_start += n_i(i);
  }

  NumericVector s_scores(n_obs);

  double a1 = sum(pow(A2,2.0));
  double a2 = -2*sum(A2*Ab);
  double a4 = sum(pow(Ab,2.0));
  double b1 = sum(A1);
  double b2 = -2*sum(Bb);
  double b4 = sum(Bc);

  double eta_optimal = 0.5*(2*a1*b2-a2*b1)/(a4*b1 - 2*a2*b2 - 2*a1*b4 + 3*(a1*pow(b2,2.0)/b1));

  return eta_optimal;
}



// Optimal step size for autoregressive proxyCCF
// [[Rcpp::export]]

double optimal_step_size_AR_cpp(NumericVector epsilon_hat, NumericVector xi_hat, int I, int n_obs, NumericVector n_i, NumericVector s, NumericVector h, double theta) {
  NumericVector A1(I);
  NumericVector A2(I);
  NumericVector Ab(I);
  NumericVector Ac(I);
  NumericVector Bb(I);
  NumericVector Bc(I);
  NumericVector alpha1(I);
  NumericVector beta1(I);
  NumericVector gamma1(I);
  NumericVector alpha2(I);
  NumericVector beta2(I);
  NumericVector gamma2(I);
  NumericVector alpha1h(I);
  NumericVector beta1h(I);
  NumericVector gamma1h(I);
  NumericVector alpha1sh(I);
  NumericVector beta1sh(I);
  NumericVector gamma1sh(I);
  NumericVector alpha2h(I);
  NumericVector beta2h(I);
  NumericVector gamma2h(I);
  NumericVector alpha2sh(I);
  NumericVector beta2sh(I);
  NumericVector gamma2sh(I);
  int n_cumsum_start = 0;
  for(int i = 0; i < I; i++) {
    NumericVector s_i(n_i(i));
    NumericVector h_i(n_i(i));
    NumericVector xi_hat_i(n_i(i));
    NumericVector epsilon_hat_i(n_i(i));
    for(int gg = 0; gg < n_i(i); gg++) {
      s_i(gg) = s(n_cumsum_start + gg);
      h_i(gg) = h(n_cumsum_start + gg);
      xi_hat_i(gg) = xi_hat(n_cumsum_start + gg);
      epsilon_hat_i(gg) = epsilon_hat(n_cumsum_start + gg);
    }
    if (n_i(i) > 2){
      NumericVector s_i_sub1(n_i(i)-1);
      NumericVector s_i_sub2(n_i(i)-1);
      NumericVector h_i_sub1(n_i(i)-1);
      NumericVector h_i_sub2(n_i(i)-1);
      NumericVector xi_hat_i_sub1(n_i(i)-1);
      NumericVector xi_hat_i_sub2(n_i(i)-1);
      NumericVector epsilon_hat_i_sub1(n_i(i)-1);
      NumericVector epsilon_hat_i_sub2(n_i(i)-1);
      for(int gg = 0; gg < (n_i(i)-1); gg++) {
        s_i_sub1(gg) = s_i(gg);
        s_i_sub2(gg) = s_i(gg+1);
        h_i_sub1(gg) = h_i(gg);
        h_i_sub2(gg) = h_i(gg+1);
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
      alpha1h(i) = sum(pow(h_i,2.0)*pow(xi_hat_i,2.0));
      beta1h(i) = alpha1(i) - pow(h_i(0)*xi_hat_i(0),2.0) - pow(h_i(n_i(i)-1)*xi_hat_i(n_i(i)-1),2.0);
      gamma1h(i) = 2*sum(h_i_sub1*h_i_sub2*xi_hat_i_sub1*xi_hat_i_sub2);
      alpha1sh(i) = sum(s_i*h_i*pow(xi_hat_i,2.0));
      beta1sh(i) = alpha1(i) - s_i(0)*h_i(0)*pow(xi_hat_i(0),2.0) - s_i(n_i(i)-1)*h_i(n_i(i)-1)*pow(xi_hat_i(n_i(i)-1),2.0);
      gamma1sh(i) = sum(h_i_sub1*s_i_sub2*xi_hat_i_sub1*xi_hat_i_sub2) + sum(s_i_sub1*h_i_sub2*xi_hat_i_sub1*xi_hat_i_sub2);
      alpha2h(i) = sum(pow(h_i,2.0)*xi_hat_i*epsilon_hat_i);
      beta2h(i) = alpha2(i) - pow(h_i(0),2.0)*xi_hat_i(0)*epsilon_hat_i(0) - pow(h_i(n_i(i)-1),2.0)*xi_hat_i(n_i(i)-1)*epsilon_hat_i(n_i(i)-1);
      gamma2h(i) = sum(h_i_sub1*h_i_sub2*(xi_hat_i_sub2*epsilon_hat_i_sub1+xi_hat_i_sub1*epsilon_hat_i_sub2));
      alpha2sh(i) = sum(s_i*h_i*xi_hat_i*epsilon_hat_i);
      beta2sh(i) = alpha2(i) - s_i(0)*h_i(0)*xi_hat_i(0)*epsilon_hat_i(0) - s_i(n_i(i)-1)*h_i(n_i(i)-1)*xi_hat_i(n_i(i)-1)*epsilon_hat_i(n_i(i)-1);
      gamma2sh(i) = sum((s_i_sub1*h_i_sub2+h_i_sub1*s_i_sub2)*0.5*(xi_hat_i_sub2*epsilon_hat_i_sub1+xi_hat_i_sub1*epsilon_hat_i_sub2));
      A1(i) = alpha1(i) + pow(theta,2.0)*beta1(i) - theta*gamma1(i);
      A2(i) = alpha2(i) + pow(theta,2.0)*beta2(i) - theta*gamma2(i);
      Ab(i) = 2*( alpha2sh(i) + pow(theta,2.0)*beta2sh(i) - theta*gamma2sh(i) );
      Ac(i) = alpha2h(i) + pow(theta,2.0)*beta2h(i) - theta*gamma2h(i);
      Bb(i) = alpha1sh(i) + pow(theta,2.0)*beta1sh(i) - theta*gamma1sh(i);
      Bc(i) = alpha1h(i) + pow(theta,2.0)*beta1h(i) - theta*gamma1h(i);
    }
    else if (n_i(i) == 2){
      alpha1(i) = sum(pow(s_i,2.0)*pow(xi_hat_i,2.0));
      gamma1(i) = 2*(s_i(1)*s_i(0)*xi_hat_i(1)*xi_hat_i(0));
      alpha2(i) = sum(pow(s_i,2.0)*xi_hat_i*epsilon_hat_i);
      gamma2(i) = s_i(1)*s_i(0)*(xi_hat_i(0)*epsilon_hat_i(1)+xi_hat_i(1)*epsilon_hat_i(0));
      alpha1h(i) = sum(pow(h_i,2.0)*pow(xi_hat_i,2.0));
      gamma1h(i) = 2*h_i(0)*h_i(1)*xi_hat_i(0)*xi_hat_i(1);
      alpha1sh(i) = sum(s_i*h_i*pow(xi_hat_i,2.0));
      gamma1sh(i) = h_i(0)*s_i(1)*xi_hat_i(0)*xi_hat_i(1) + s_i(0)*h_i(1)*xi_hat_i(0)*xi_hat_i(1);
      alpha2h(i) = sum(pow(h_i,2.0)*xi_hat_i*epsilon_hat_i);
      gamma2h(i) = h_i(0)*h_i(1)*(xi_hat_i(1)*epsilon_hat_i(0)+xi_hat_i(0)*epsilon_hat_i(1));
      alpha2sh(i) = sum(s_i*h_i*xi_hat_i*epsilon_hat_i);
      gamma2sh(i) = (s_i(0)*h_i(1)+h_i(0)*s_i(1))*0.5*(xi_hat_i(1)*epsilon_hat_i(0)+xi_hat_i(0)*epsilon_hat_i(1));
      A1(i) = alpha1(i) - theta*gamma1(i);
      A2(i) = alpha2(i) - theta*gamma2(i);
      Ab(i) = 2*( alpha2sh(i) - theta*gamma2sh(i) );
      Ac(i) = alpha2h(i) - theta*gamma2h(i);
      Bb(i) = alpha1sh(i) - theta*gamma1sh(i);
      Bc(i) = alpha1h(i) - theta*gamma1h(i);
    }
    else if (n_i(i) == 1){
      alpha1(i) = 2*pow(s_i(0),2.0)*pow(xi_hat_i(0),2.0);
      alpha2(i) = 2*pow(s_i(0),2.0)*xi_hat_i(0)*epsilon_hat_i(0);
      alpha1h(i) = pow(h_i(0),2.0)*pow(xi_hat_i(0),2.0);
      alpha1sh(i) = s_i(0)*h_i(0)*pow(xi_hat_i(0),2.0);
      alpha2h(i) = pow(h_i(0),2.0)*xi_hat_i(0)*epsilon_hat_i(0);
      alpha2sh(i) = s_i(0)*h_i(0)*xi_hat_i(0)*epsilon_hat_i(0);
      A1(i) = (1-pow(theta,2.0))*alpha1(i);
      A2(i) = (1-pow(theta,2.0))*alpha2(i);
      Ab(i) = (1-pow(theta,2.0))*2*alpha2sh(i);
      Ac(i) = (1-pow(theta,2.0))*alpha2h(i);
      Bb(i) = (1-pow(theta,2.0))*alpha1sh(i);
      Bc(i) = (1-pow(theta,2.0))*alpha1h(i);
    }
    n_cumsum_start += n_i(i);
  }

  double a1 = sum(pow(A2,2.0));
  double a2 = -2*sum(A2*Ab);
  double a4 = sum(pow(Ab,2.0));
  double b1 = sum(A1);
  double b2 = -2*sum(Bb);
  double b4 = sum(Bc);

  double eta_optimal = 0.5*(2*a1*b2-a2*b1)/(a4*b1 - 2*a2*b2 - 2*a1*b4 + 3*(a1*pow(b2,2.0)/b1));

  return eta_optimal;
}










