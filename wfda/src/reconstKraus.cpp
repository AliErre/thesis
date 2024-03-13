#include "reconstKraus.h"
//reconstructs the curve at index
List reconstKraus_fun(const NumericMatrix& Y, const NumericVector& mean_vec,
                      const NumericMatrix& cov_mat,unsigned index, double alpha) {
  //index = reconst_fcts[i], o in gcv è un'altra cosa

  arma::mat cov_mat_arma = as<arma::mat>(cov_mat);
  NumericVector X_cent_vec = Y(_,index) - mean_vec;//Y = m_Y, mean_vec = m_mean
  arma::vec X_cent_vec_arma = as<arma::vec>(X_cent_vec);//as<arma::vec> converts NA_REAL to "Nan"
  arma::uvec M_bool = arma::find_nonfinite(X_cent_vec_arma);//finds NaN
  arma::uvec O_bool = arma::find_finite(X_cent_vec_arma); 

  arma::mat covMO_mat = cov_mat_arma.submat(M_bool, O_bool);  //cov_mat = m_cov
  arma::mat covOO_mat = cov_mat_arma.submat(O_bool, O_bool);

  arma::mat covOO_a_mat = covOO_mat + alpha * arma::eye<arma::mat>(O_bool.n_elem, O_bool.n_elem);
  arma::mat covOO_a_mat_inv = arma::inv(covOO_a_mat);
  arma::mat Aa_mat = covMO_mat * covOO_a_mat_inv;
  arma::vec X_M_fit_cent_vec = Aa_mat * X_cent_vec_arma.rows(O_bool);

  arma::vec X_cent_reconst_vec = X_cent_vec_arma;
  X_cent_reconst_vec.rows(M_bool) = X_M_fit_cent_vec;
  arma::vec eigenvalues = arma::eig_sym(covOO_mat);
  arma::vec lam00 = eigenvalues.elem(arma::find(eigenvalues > 0));
  double df = arma::sum(lam00/(lam00 + alpha));

  //compute h(t)
  arma::mat covMM_mat = cov_mat_arma.submat(M_bool, M_bool);
  arma::mat V_mat = cov_MM_mat - Aa_mat * covOO_mat * Aa_mat.t();

  arma::vec hi;
  if(all(arma::diagvec(V_mat) > 0))
  {
    arma::vec vi = arma::sqrt(arma::diagvec(V_mat));
    double h0 = 0.2*max(vi);

    hi = vi;
    for (uword i = 0; i < hi.n_elem; ++i) {
      hi(i) = h0 < vi(i) ? vi(i) : h0;
    }}
  else{
    hi = arma::vec(V_mat.n_rows, fill::ones) * 0.01;
  }

  arma::vec hi_scaled;
  if(hi.n_rows == 1){
    hi_scaled = hi/arma::sqrt(covMM_mat);
  }
  else{
    hi_scaled = hi/arma::sqrt(arma::diagvec(covMM_mat));   
  }
  return List::create(_["X_cent_reconst_vec"] = X_cent_reconst_vec, //arma::vec
                      _["df"] = df,//double
                      _["alpha"] = alpha,//double
                      _["hi"] = hi_scaled);//arma::vec
}