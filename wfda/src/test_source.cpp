#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
//[[Rcpp::export]]
arma::uvec test(const NumericVector& M)
{
    arma::uvec v = arma::find_nonfinite(as<arma::vec>(M));
    return v;
}