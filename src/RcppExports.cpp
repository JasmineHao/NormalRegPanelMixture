// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cppRegPanelmixPMLE
List cppRegPanelmixPMLE(NumericMatrix bs, NumericVector ys, NumericMatrix xs, NumericMatrix zs, NumericVector mu0s, NumericVector sigma0s, int m, int p, int t, double an, int maxit, int ninits, double tol, double tau, int h, int k);
RcppExport SEXP _NormalRegPanelMixture_cppRegPanelmixPMLE(SEXP bsSEXP, SEXP ysSEXP, SEXP xsSEXP, SEXP zsSEXP, SEXP mu0sSEXP, SEXP sigma0sSEXP, SEXP mSEXP, SEXP pSEXP, SEXP tSEXP, SEXP anSEXP, SEXP maxitSEXP, SEXP ninitsSEXP, SEXP tolSEXP, SEXP tauSEXP, SEXP hSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type bs(bsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ys(ysSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type xs(xsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type zs(zsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu0s(mu0sSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma0s(sigma0sSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type an(anSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int >::type ninits(ninitsSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(cppRegPanelmixPMLE(bs, ys, xs, zs, mu0s, sigma0s, m, p, t, an, maxit, ninits, tol, tau, h, k));
    return rcpp_result_gen;
END_RCPP
}
// cppnormalpanelmixPMLE
List cppnormalpanelmixPMLE(NumericMatrix bs, NumericVector ys, NumericMatrix zs, NumericVector mu0s, NumericVector sigma0s, int m, int p, int t, double an, int maxit, int ninits, double tol, double tau, int h, int k);
RcppExport SEXP _NormalRegPanelMixture_cppnormalpanelmixPMLE(SEXP bsSEXP, SEXP ysSEXP, SEXP zsSEXP, SEXP mu0sSEXP, SEXP sigma0sSEXP, SEXP mSEXP, SEXP pSEXP, SEXP tSEXP, SEXP anSEXP, SEXP maxitSEXP, SEXP ninitsSEXP, SEXP tolSEXP, SEXP tauSEXP, SEXP hSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type bs(bsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ys(ysSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type zs(zsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu0s(mu0sSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma0s(sigma0sSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type an(anSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int >::type ninits(ninitsSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(cppnormalpanelmixPMLE(bs, ys, zs, mu0s, sigma0s, m, p, t, an, maxit, ninits, tol, tau, h, k));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NormalRegPanelMixture_cppRegPanelmixPMLE", (DL_FUNC) &_NormalRegPanelMixture_cppRegPanelmixPMLE, 16},
    {"_NormalRegPanelMixture_cppnormalpanelmixPMLE", (DL_FUNC) &_NormalRegPanelMixture_cppnormalpanelmixPMLE, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_NormalRegPanelMixture(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
