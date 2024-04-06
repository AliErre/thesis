#ifndef RECONST_KRAUS_H
#define RECONST_KRAUS_H
#include <RcppArmadillo.h>
using namespace Rcpp;
List reconstKraus_fun(const NumericMatrix&, const NumericVector&,
                      const NumericMatrix&, unsigned, double);

#endif