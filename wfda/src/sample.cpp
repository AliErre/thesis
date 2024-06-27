#include <RcppArmadilloExtensions/sample.h>
arma::vec sample(const arma::vec& x, const size_t& size, const bool& replace){
  return Rcpp::RcppArmadillo::sample(x, size, replace);//forse int anziche size_t
}
arma::uvec sample(const arma::uvec& x, const size_t& size, const bool& replace){
  return Rcpp::RcppArmadillo::sample(x, size, replace);
}