#ifndef RECONST_KRAUS_H
#define RECONST_KRAUS_H
#include <Rcpp.h>
#include<RcppArmadillo.h>
List reconstKraus_fun(const NumericMatrix& Y, const NumericVector& mean_vec,
                      const NumericMatrix& cov_mat,unsigned index, double alpha){};

#endif RECONST_KRAUS_H