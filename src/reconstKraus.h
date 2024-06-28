#ifndef RECONST_KRAUS_H
#define RECONST_KRAUS_H
#include <RcppArmadillo.h>
using namespace Rcpp;
std::tuple<NumericVector, double, double, arma::vec, arma::uvec> 
reconstKraus_fun(const NumericMatrix&, const NumericVector&, const NumericMatrix&, unsigned, double);

#endif