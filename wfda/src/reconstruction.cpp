#include "reconstruction.h"
#include <algorithm>

std::vector<int> find_obs_inc(const NumericMatrix& Y) const{//forse questo non metterlo come membro della classe ma FREE FUNCTION
//se lo metto come free function potrei chiamarlo a prescindere dall'avere un oggetto della classe
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


const std::vector<double>& meanKraus() {

  int nRows = m_Y.row();
  int nCols = m_Y.col();
  m_mean.resize(nRows);
    for (int i = 0; i < nRows; ++i) {
        double sum = 0;
        int naCount = 0;

        NumericMatrix::ConstRow row = X_mat(i, _); //::COnstRow gives constant reference to the current row
        // Iterate over the elements of the row
        for (auto it = row.begin(); it != row.end(); ++it) {
            if (NumericVector::is_na(*it)) {
                naCount++;
            } else {
                sum += *it;
            }
        }
        
        m_mean[i] = naCount < nCols ? sum / (nCols - naCount) : NA_REAL;
    }

    return m_mean;
}

const NumericMatrix& covKraus(){
  int nRows = m_Y.row();
  int nCols = m_Y.col();

  std::vector<double> rowMeans = m_mean.empty() ? meanKraus(m_Y) : m_mean;
  NumericMatrix X_cent_mat(nRows, nCols);//fixed dimensions 

    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            X_cent_mat(i, j) = NumericVector::is_na(m_Y(i, j)) ? NA_REAL : m_Y(i, j) - rowMeans[i];//ricontrolla
        }
    }
    
    //account for NAs
    m_cov = NumericMatrix(nRows, nRows);
    for (int s = 0; s < nRows; ++s) {
        for (int t = s; t < nRows; ++t) {
            double sum = 0.0;
            int count = 0;
            for (int k = 0; k < nCols; ++k) {
                if (!NumericVector::is_na(X_cent_mat(s, k)) && !NumericVector::is_na(X_cent_mat(t, k))) {
                    sum += X_cent_mat(s, k) * X_cent_mat(t, k);
                    ++count;
                }
            }
            double covValue = count > 0 ? sum / count : NA_REAL;
            m_cov(s, t) = covValue;
            m_cov(t, s) = covValue; //symmetric
        }
    }
    
  return m_cov;
}

double gcvKraus(const std::vector<std::vector<double>>& covMat, const std::vector<double>& meanVec, 
                        const NumericMatrix& X, const bool M_bool_vec, const Numeric& alpha) const{

                        }