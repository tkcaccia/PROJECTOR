// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393
arma::vec SVDS(arma::sp_mat  &X)


#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


// SVDS
arma::vec SVDS(arma::sp_mat  &X);
RcppExport SEXP PROJECTOR(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::sp_mat >::type X(XSEXP);
    __result = Rcpp::wrap(SVDS(X));
    return __result;
END_RCPP
}
