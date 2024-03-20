#include "reconstruction.h"
#include <algorithm>
#include <numeric>
#include <limits>
#include "gcv.h"

//capire se mi serve la forward declaration di reconstKraus
//destructor per gcv?

std::vector<int> ReconstructionBase::find_obs_inc(const NumericMatrix& Y) const{
//forse questo non metterlo come membro della classe ma FREE FUNCTION
//se lo metto come free function potrei chiamarlo a prescindere dall'avere un oggetto della classe
  int n = Y.ncol(); 
  std::vector<int> reconst_fcts; //initialized empty. mettere unsigned?
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


IntegerVector ReconstructionBase::reconst_fcts() const{
  //when return to R, indices must start from 1
  std::vector<int> shifted_reconst_fcts(m_reconst_fcts.size());
  std::transform(m_reconst_fcts.begin(), m_reconst_fcts.end(), shifted_reconst_fcts.begin(),
                 [](int value){return value + 1; }); //lambda function
  return wrap(shifted_reconst_fcts); //wrap does not perform additional memory allocation, directly references
}


const NumericVector& ReconstructionBase::meanRows() {

  int nRows = m_Y.nrow();
  int nCols = m_Y.ncol();
  m_mean = NumericVector(nRows);//vedi se lasciarlo o no
    for (int i = 0; i < nRows; ++i) {
        double sum = 0;
        int naCount = 0;

        NumericVector row = m_Y(i, _); //::ConstRow gives constant reference to the current row
        // Iterate over the elements of the row
        for (auto it = row.begin(); it != row.end(); ++it) {
            if (NumericVector::is_na(*it)) {
                naCount++;
            } else {
                sum += *it;
            }
        }
        m_mean[i] = sum / (nCols - naCount);
        /*m_mean[i] = naCount < nCols ? sum / (nCols - naCount) : NA_REAL; 
        //in realtà non dovrebbe mai capitare di avere degli NA
        //ci dovrebbe essere sempre almeno una curva completa*/
    }

    return m_mean;//retunr a reference to the data member
}

const NumericMatrix& ReconstructionBase::covMatrix(){
  //if(m_mean.length() == 0){meanRows();} not necessary
  int nRows = m_Y.nrow();
  int nCols = m_Y.ncol();
  NumericMatrix X_cent_mat(nRows, nCols);//fixed dimensions 

    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            X_cent_mat(i, j) = NumericVector::is_na(m_Y(i, j)) ? NA_REAL : (m_Y(i, j) - m_mean[i]);
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
            //double covValue = count > 0 ? sum / count : NA_REAL; commented out perchè penso che non possa mai capitare di avere NA nella cov_mat
            double covValue = sum/count;
            m_cov(s, t) = covValue;
            m_cov(t, s) = covValue; //symmetric
        }
    }
    
  return m_cov;
}


List ReconstructionKraus::reconstructCurve(double alpha = 0.0, bool all = FALSE,const NumericVector& periods = NumericVector()) const{
//dovrei avere già mean_vec e cov_mat nella classe appena chiamo il costruttore
  int n = m_Y.ncol();
  IntegerVector reconst_fcts;
  //if all = TRUE -> ricostruiscile tutte
  if(all){
    reconst_fcts = seq_len(n); //Rcpp sugar, vedi se l'alternativa con for è più efficiente
  }
  else{
    reconst_fcts = (*this).reconst_fcts() - 1;//getter, trasformava gli indici avanti di uno perchè li restituiva ad R. Quindi ora togli 1
  }
  int r = m_Y.nrow();
  NumericMatrix X_reconst_mat(r,reconst_fcts.size());
  int column = 0;
  for(int& index: reconst_fcts){
    X_reconst_mat(_,column) = m_Y(_,index);
    column++;
  }
  NumericMatrix W_reconst_mat(r,reconst_fcts.length());//initialized filled with 0s
  //std::fill(W_reconst_mat.begin(),W_reconst_mat.end(),1);//posso fare così grazie a come sono salvate le NumericMatrix in memoria

  std::vector<int> nonNA_fcts; //mask in R
  nonNA_fcts.reserve(n);//metterlo come data member?
  //apply(X_mat,2,function(x)!any(is.na(x)))
  for(int i = 0; i < n; i++){
    //trasforma tutti i loop con size_t
    NumericVector col = m_Y(_, i);
    if(is_false(any(is_na(col))))
      nonNA_fcts.push_back(i); 
  }

  NumericMatrix X_Compl_mat(r,nonNA_fcts.size());
  column = 0;
  for(auto& index: nonNA_fcts){
    X_Compl_mat(_,column) = m_Y(_,index);
    column++;
  }

  NumericVector alpha_vec(reconst_fcts.size());
  NumericVector df_vec(reconst_fcts.size());
  
  column = 0;
  gcv GCV(X_Compl_mat, m_mean, m_cov);//vedere se poi chiamare destructor. Only created once
  
  for(auto& index:reconst_fcts){

    LogicalVector M_bool = is_na(m_Y(_,index));
    LogicalVector O_bool = !M_bool;
    GCV.set_bool(M_bool);
    
    if(alpha == 0.0)//R_NilValue = NULL in R, in attesa di capire come gestire NULL, metto 0.0 di default
    {
      double max_bound = 0.0;
      for(int i = 0; i < r;++i){//r should be the nrow, ncol of m_cov -> diag is of length r
        max_bound += m_cov(i,i);
      }
      alpha_vec[column] = optimize(&GCV, std::numeric_limits<double>::epsilon(), max_bound, false, 1e-3); //false -> minimization
    }else{
      alpha_vec[column] = alpha;
    }

    List resultKraus = reconstKraus_fun(m_Y, m_mean, m_cov,index, alpha_vec[column]);//forse dovrei cambiarla e passargli M_bool
  //reconstKraus["X_cent_reconst_vec"] is an arma::vec
    NumericVector X_reconst = resultKraus["X_cent_reconst_vec"];
    X_reconst_mat(_,column) = X_reconst + m_mean;
    df_vec.push_back(resultKraus["df"]);//check
    // W_reconst_mat[M_bool_vec,i]
    NumericVector hi = resultKraus["hi"];
    for(int i = 0; i < r;i++){
      if(M_bool[i])//gli altri pesi rimangono ad 1
        W_reconst_mat(i,column) = 1 - hi[i-hi.length()-1];//controlla hi[] che sia giusto
      else
        W_reconst_mat(i,column) = 1;
    }
    column++;
  }
  return(List::create(_["X_reconst_mat"] = X_reconst_mat,
                      _["alpha"]         = alpha_vec, 
                      _["df"]            = df_vec,
                      _["W_reconst_mat"] = W_reconst_mat)); 

}


//extrapolation method for reconstruction from last observed period
List ReconstructionExtrapolation::reconstructCurve(double alpha = 0.0, bool all = FALSE,const NumericVector& periods = NumericVector()) const{
  int r = periods.length(); //m_Y.nrow()
  int n = m_Y.ncol();
  double sum = 0.0;
  NumericVector mean_slope(r-1); //slope media per ogni riga
  IntegerVector no_rec_fcts(n-m_reconst_fcts.size());
  //build vector of indeces not in m_reconst_fcts
  int row = 0;
  double sum = 0.0;
  for( int i = 0; i < m_Y.ncol();i++){
    bool found = std::binary_search(m_reconst_fcts.begin(),m_reconst_fcts.end(),i);//works cause reconst_fcts is an ordered vector
    if(!found){no_rec_fcts(row) = i; row++;}
  }

  for(int j = 0; j < r-1; j++)
  {
    for(auto& index:no_rec_fcts)
      sum += m_Y(r-1,index) - m_Y(j,index);
    mean_slope(j) = sum/(no_rec_fcts.length()*(periods(r) - periods(j)));
  }

  NumericMatrix Y_reconstruct(r,n);
  for(auto& index:m_reconst_fcts)
  {
    const NumericVector& col_Y = m_Y(_,index);//trasforma tutte le copie che hai fatto in const references per EVITARE COPIE
    const LogicalVector& id_na = is_na(col_Y);//check nelle slide pacs se posso bindare const ref a un temporary
    NumericVector observed_periods = periods[!id_na];
    double last_obs_period = observed_periods[observed_periods.size()-1];
    const NumericVector& slopes = mean_slope[!id_na];
    double last_slope = slopes[slopes.size()-1];
    const NumericVector& observed_col_Y = col_Y[!id_na];
    double last_obs_value = observed_col_Y[observed_col_Y.size()-1];
    //reconstruct with equation of a straight line
    for(auto& index_na:id_na){
      if(index_na)
        Y_reconstruct(index_na, index) = last_obs_value + last_slope*(periods(index_na)-last_obs_period);
    }
  }
  return List::create(_["Extrapolated_curves"] = Y_reconstruct,
                      _["mean_slopes"] = mean_slope);

}
