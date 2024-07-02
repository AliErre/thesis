#ifndef SAMPLE_H
#define SAMPLE_H

#include <RcppArmadillo.h>
arma::vec sample(const arma::vec& x, const size_t& size, const bool& replace);
arma::uvec sample(const arma::uvec& x, const size_t& size, const bool& replace);
#endif