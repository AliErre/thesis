#ifndef KNEIP_LIEBL_H
#define KNEIP_LIEBL_H
#include <RcppArmadillo.h>
#include <vector>
using namespace Rcpp;
//inner reconstruction function called in the gcv, with gcv parameter k
//if method is KLAl5, fragmO.size() = 0
List reconstKL_fun(const NumericVector&, const std::vector<double>&, const arma::uvec&, 
                   const arma::vec&, const NumericMatrix&, const NumericVector&, int);

std::pair<std::vector<double>,NumericMatrix> irreg2mat(const std::vector<std::tuple<int, double, double>>&, 
                                                       bool binning, int max_bins);

std::pair<NumericMatrix,NumericVector> smooth_cov(const NumericMatrix&,const NumericMatrix&,const NumericVector&,
                                                  int, int, int);

std::vector<double> quadWeights(const std::vector<double>&, const std::string& method = "trapezoidal");

double weighted_mean(const std::vector<double>&, const std::vector<double>&);

double trapezioidal_rule(const arma::vec&, const arma::vec&);

std::tuple<List, List, List, arma::mat, List, List, arma::vec, List, List, List, double, arma::mat> 
eigen(const std::vector<double>&, const std::vector<std::vector<double>>&,const NumericMatrix&, 
      double, const NumericMatrix&, const NumericVector&, const NumericVector&, 
      bool, const IntegerVector&);

int gcvKneipLiebl(const NumericVector&, const std::pair<std::vector<double>, NumericMatrix>&, 
                  const std::vector<double>&, const NumericVector&, const arma::uvec&,
                  const arma::mat&, double, const std::string&, double pev = 0.99);
#endif



