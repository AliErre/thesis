#include "reconstruction.h"

std::vector<size_t> find_obs_inc(const NumericMatrix& Y) {
  int n = Y.ncol(); 
  std::vector<size_t> reconst_fcts;
  reconst_fcts.reserve(n); //reserve memory -> makes push_back better

  //syntax that follows is allowed thanks to Rcpp
  for (int i = 0; i < n; ++i) {
    NumericVector col = Y(_, i);
    // Check NA
    if (is_true(any(is_na(col)))) {
      reconst_fcts.push_back(i); //0-based index
    }
  }
  
  reconst_fcts.shrink_to_fit(); //release memory

  return reconst_fcts;
}

