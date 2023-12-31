// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// optimal_step_size_equicorr_cpp
double optimal_step_size_equicorr_cpp(NumericVector epsilon_hat, NumericVector xi_hat, int I, int n_obs, NumericVector n_i, NumericVector s, NumericVector h, double theta);
RcppExport SEXP _sandwich_boost_optimal_step_size_equicorr_cpp(SEXP epsilon_hatSEXP, SEXP xi_hatSEXP, SEXP ISEXP, SEXP n_obsSEXP, SEXP n_iSEXP, SEXP sSEXP, SEXP hSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type epsilon_hat(epsilon_hatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xi_hat(xi_hatSEXP);
    Rcpp::traits::input_parameter< int >::type I(ISEXP);
    Rcpp::traits::input_parameter< int >::type n_obs(n_obsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_i(n_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(optimal_step_size_equicorr_cpp(epsilon_hat, xi_hat, I, n_obs, n_i, s, h, theta));
    return rcpp_result_gen;
END_RCPP
}
// optimal_step_size_AR_cpp
double optimal_step_size_AR_cpp(NumericVector epsilon_hat, NumericVector xi_hat, int I, int n_obs, NumericVector n_i, NumericVector s, NumericVector h, double theta);
RcppExport SEXP _sandwich_boost_optimal_step_size_AR_cpp(SEXP epsilon_hatSEXP, SEXP xi_hatSEXP, SEXP ISEXP, SEXP n_obsSEXP, SEXP n_iSEXP, SEXP sSEXP, SEXP hSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type epsilon_hat(epsilon_hatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xi_hat(xi_hatSEXP);
    Rcpp::traits::input_parameter< int >::type I(ISEXP);
    Rcpp::traits::input_parameter< int >::type n_obs(n_obsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_i(n_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(optimal_step_size_AR_cpp(epsilon_hat, xi_hat, I, n_obs, n_i, s, h, theta));
    return rcpp_result_gen;
END_RCPP
}
// risk_equicorr_cpp
double risk_equicorr_cpp(NumericVector epsilon_hat, NumericVector xi_hat, int I, int n_obs, NumericVector n_i, NumericVector s, double theta);
RcppExport SEXP _sandwich_boost_risk_equicorr_cpp(SEXP epsilon_hatSEXP, SEXP xi_hatSEXP, SEXP ISEXP, SEXP n_obsSEXP, SEXP n_iSEXP, SEXP sSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type epsilon_hat(epsilon_hatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xi_hat(xi_hatSEXP);
    Rcpp::traits::input_parameter< int >::type I(ISEXP);
    Rcpp::traits::input_parameter< int >::type n_obs(n_obsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_i(n_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(risk_equicorr_cpp(epsilon_hat, xi_hat, I, n_obs, n_i, s, theta));
    return rcpp_result_gen;
END_RCPP
}
// risk_AR_cpp
double risk_AR_cpp(NumericVector epsilon_hat, NumericVector xi_hat, int I, int n_obs, NumericVector n_i, NumericVector s, double theta);
RcppExport SEXP _sandwich_boost_risk_AR_cpp(SEXP epsilon_hatSEXP, SEXP xi_hatSEXP, SEXP ISEXP, SEXP n_obsSEXP, SEXP n_iSEXP, SEXP sSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type epsilon_hat(epsilon_hatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xi_hat(xi_hatSEXP);
    Rcpp::traits::input_parameter< int >::type I(ISEXP);
    Rcpp::traits::input_parameter< int >::type n_obs(n_obsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_i(n_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(risk_AR_cpp(epsilon_hat, xi_hat, I, n_obs, n_i, s, theta));
    return rcpp_result_gen;
END_RCPP
}
// risk_nested_cpp
double risk_nested_cpp(NumericVector epsilon_hat, NumericVector xi_hat, NumericVector J, int n_obs, NumericMatrix n_ij, NumericVector s, NumericVector theta);
RcppExport SEXP _sandwich_boost_risk_nested_cpp(SEXP epsilon_hatSEXP, SEXP xi_hatSEXP, SEXP JSEXP, SEXP n_obsSEXP, SEXP n_ijSEXP, SEXP sSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type epsilon_hat(epsilon_hatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xi_hat(xi_hatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type J(JSEXP);
    Rcpp::traits::input_parameter< int >::type n_obs(n_obsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type n_ij(n_ijSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(risk_nested_cpp(epsilon_hat, xi_hat, J, n_obs, n_ij, s, theta));
    return rcpp_result_gen;
END_RCPP
}
// ngradient_equicorr_cpp
List ngradient_equicorr_cpp(NumericVector epsilon_hat, NumericVector xi_hat, int I, int n_obs, NumericVector n_i, NumericVector s, double theta);
RcppExport SEXP _sandwich_boost_ngradient_equicorr_cpp(SEXP epsilon_hatSEXP, SEXP xi_hatSEXP, SEXP ISEXP, SEXP n_obsSEXP, SEXP n_iSEXP, SEXP sSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type epsilon_hat(epsilon_hatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xi_hat(xi_hatSEXP);
    Rcpp::traits::input_parameter< int >::type I(ISEXP);
    Rcpp::traits::input_parameter< int >::type n_obs(n_obsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_i(n_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(ngradient_equicorr_cpp(epsilon_hat, xi_hat, I, n_obs, n_i, s, theta));
    return rcpp_result_gen;
END_RCPP
}
// ngradient_AR_cpp
List ngradient_AR_cpp(NumericVector epsilon_hat, NumericVector xi_hat, int I, int n_obs, NumericVector n_i, NumericVector s, double theta);
RcppExport SEXP _sandwich_boost_ngradient_AR_cpp(SEXP epsilon_hatSEXP, SEXP xi_hatSEXP, SEXP ISEXP, SEXP n_obsSEXP, SEXP n_iSEXP, SEXP sSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type epsilon_hat(epsilon_hatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xi_hat(xi_hatSEXP);
    Rcpp::traits::input_parameter< int >::type I(ISEXP);
    Rcpp::traits::input_parameter< int >::type n_obs(n_obsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_i(n_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(ngradient_AR_cpp(epsilon_hat, xi_hat, I, n_obs, n_i, s, theta));
    return rcpp_result_gen;
END_RCPP
}
// ngradient_nested_cpp
List ngradient_nested_cpp(NumericVector epsilon_hat, NumericVector xi_hat, NumericVector J, int n_obs, NumericMatrix n_ij, NumericVector s, NumericVector theta);
RcppExport SEXP _sandwich_boost_ngradient_nested_cpp(SEXP epsilon_hatSEXP, SEXP xi_hatSEXP, SEXP JSEXP, SEXP n_obsSEXP, SEXP n_ijSEXP, SEXP sSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type epsilon_hat(epsilon_hatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xi_hat(xi_hatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type J(JSEXP);
    Rcpp::traits::input_parameter< int >::type n_obs(n_obsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type n_ij(n_ijSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(ngradient_nested_cpp(epsilon_hat, xi_hat, J, n_obs, n_ij, s, theta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sandwich_boost_optimal_step_size_equicorr_cpp", (DL_FUNC) &_sandwich_boost_optimal_step_size_equicorr_cpp, 8},
    {"_sandwich_boost_optimal_step_size_AR_cpp", (DL_FUNC) &_sandwich_boost_optimal_step_size_AR_cpp, 8},
    {"_sandwich_boost_risk_equicorr_cpp", (DL_FUNC) &_sandwich_boost_risk_equicorr_cpp, 7},
    {"_sandwich_boost_risk_AR_cpp", (DL_FUNC) &_sandwich_boost_risk_AR_cpp, 7},
    {"_sandwich_boost_risk_nested_cpp", (DL_FUNC) &_sandwich_boost_risk_nested_cpp, 7},
    {"_sandwich_boost_ngradient_equicorr_cpp", (DL_FUNC) &_sandwich_boost_ngradient_equicorr_cpp, 7},
    {"_sandwich_boost_ngradient_AR_cpp", (DL_FUNC) &_sandwich_boost_ngradient_AR_cpp, 7},
    {"_sandwich_boost_ngradient_nested_cpp", (DL_FUNC) &_sandwich_boost_ngradient_nested_cpp, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_sandwich_boost(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
