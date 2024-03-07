#include "reconstruction.h"
#include <algorithm>

std::vector<int> find_obs_inc(const NumericMatrix& Y) const{
  int n = Y.ncol(); 
  std::vector<int> reconst_fcts; //initialized empty
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
  return reconst_fcts; //empty if no NAs found
}


IntegerVector reconst_fcts() const{
  //when return to R, indices must start from 1
  std::vector<int> shifted_reconst_fcts(m_reconst_fcts.size());
  std::transform(m_reconst_fcts.begin(), m_reconst_fcts.end(), shifted_reconst_fcts.begin(),
                 [](int value){return value + 1; }); //lambda function
  return wrap(shifted_reconst_fcts); //wrap does not perform additional memory allocation, directly references
}


const std::vector<double>& meanKraus(const NumericMatrix& X_mat) {//oppure farla void?
  int nRows = X_mat.nrow();
  int nCols = X_mat.ncol();
  std::vector<double> rowMeans(nRows);

  for (int i = 0; i < nRows; ++i) {
    double sum = 0;
    int naCount = 0;
    // Obtain iterators to the beginning and end of the current row
    NumericMatrix::ConstRow row = X_mat(i, _);
    for (auto it = row.begin(); it != row.end(); ++it) {
      if (NumericVector::is_na(*it)) {
        naCount++;
      } else {
        sum += *it;
      }
    }
    // Calculate mean if there are non-NA values, else NA
    rowMeans[i] = naCount < nCols ? sum / (nCols - naCount) : NA_REAL;
  }
  m_mean = rowMeans
  return m_mean
}

//getter
NumericMatrix cov() const {//since there's no wrap() to go from std::vector<std::vector<double>> to NumericMatrix
    NumericMatrix Cov(m_cov.size(), m_cov[0].size());
    for (int i = 0; i < m_cov.size(); ++i) {
        for (int j = 0; j < m_cov[0].size(); ++j) {
            Cov(i, j) = m_cov[i][j];
        }
    }
    return Cov; 
}