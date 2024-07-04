#include "reconstruction.h"
#include <algorithm>
#include <numeric>
#include <limits>
#include <list>
#include <set>
#include <unordered_set>
#include <utility>
#include <stdexcept>
#include "gcv.h"
#include "KneipLiebl.h"
#include "helper_functions.h"
#include <omp.h>

//destructor per gcv?
//oppure ritorna un pair<vettore di classi, matrice Y>, vettore di classi sarrebbe colname ed è lower bound di classes

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
  //getter to be called from R: indices must start from 1
  IntegerVector rf(m_reconst_fcts.begin(), m_reconst_fcts.end());
  return rf + 1;
}


const NumericVector& ReconstructionBase::meanRows() {

  int nRows = m_Y.nrow();
  int nCols = m_Y.ncol();
  m_mean = NumericVector(nRows);
  #pragma omp parallel for if(nRows > 20)    
    for (int i = 0; i < nRows; ++i) {
        double sum = 0;
        int naCount = 0;

        NumericVector row = m_Y(i, _); 
        for (auto it = row.begin(); it != row.end(); ++it) {
            if (NumericVector::is_na(*it)) {
                naCount++;
            } else {
                sum += *it;
            }
        }
        m_mean[i] = sum / (nCols - naCount);
    }

    return m_mean;
}

const NumericMatrix& ReconstructionBase::covMatrix(){
  //if(m_mean.length() == 0){meanRows();} not necessary
  int nRows = m_Y.nrow();
  int nCols = m_Y.ncol();
  NumericMatrix X_cent_mat(nRows, nCols); 
  #pragma omp parallel for collapse(2)
    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            X_cent_mat(i, j) = NumericVector::is_na(m_Y(i, j)) ? 
                             NA_REAL : (m_Y(i, j) - m_mean[i]);
        }
    }
    
    m_cov = NumericMatrix(nRows, nRows);
    for (int s = 0; s < nRows; ++s) {
      for (int t = s; t < nRows; ++t) {
        double sum = 0.0;
        int count = 0;

        #pragma omp parallel for reduction(+:sum, count)
        for (int k = 0; k < nCols; ++k) {
          if (!NumericVector::is_na(X_cent_mat(s, k)) && !NumericVector::is_na(X_cent_mat(t, k))) {
             sum += X_cent_mat(s, k) * X_cent_mat(t, k);
              ++count;
            }
          }   
        double covValue = sum/count;
        m_cov(s, t) = covValue;
        m_cov(t, s) = covValue; //symmetric
        }     
    }
    
  return m_cov;
}



List ReconstructionKraus::reconstructCurve(Nullable<double> alpha_nullable = R_NilValue, bool all = FALSE, const Nullable<NumericVector>& t_points = R_NilValue, Nullable<int> K = R_NilValue, Nullable<int> maxBins = R_NilValue, Nullable<int> nRegGrid = R_NilValue) {

  double alpha;
  if(alpha_nullable.isNotNull()){alpha = as<double>(alpha_nullable);} 
  size_t n = m_Y.ncol();
  size_t r = m_Y.nrow();
  IntegerVector reconst_fcts;
  
  if(all){//reconstruct all
    reconst_fcts = seq_len(n) - 1; //Rcpp sugar
  }
  else{//only reconstruct the partially observed
    reconst_fcts = (*this).reconst_fcts() - 1;//getter
  }

  m_length_reconst_fcts = reconst_fcts.size();

  NumericMatrix X_reconst_mat(r,m_length_reconst_fcts);
  int column = 0;
  for(int& index: reconst_fcts){
    X_reconst_mat(_,column++) = m_Y(_,index);
  }

  NumericMatrix W_reconst_mat(r,reconst_fcts.length());//initialized with 0s
  std::fill(W_reconst_mat.begin(),W_reconst_mat.end(),1);//thanks to how NumericMatrix stored in memory

  std::vector<size_t> nonNA_fcts; //mask in R
  nonNA_fcts.reserve(n);

  for(size_t i = 0; i < n; i++){
    //trasforma tutti i loop con size_t
    NumericVector col = m_Y(_, i);
    if(is_false(any(is_na(col))))
      nonNA_fcts.push_back(i); 
  }

  NumericMatrix X_Compl_mat(r,nonNA_fcts.size());
  column = 0;
  for(const auto& index: nonNA_fcts){
    X_Compl_mat(_,column++) = m_Y(_,index);
  }

  NumericVector alpha_vec(m_length_reconst_fcts);//chose Rcpp vectors cause they have to be returned. No overhead
  NumericVector df_vec(m_length_reconst_fcts);
  
  column = 0;
  gcv GCV(X_Compl_mat, m_mean, m_cov);//vedere se poi chiamare destructor. Only created once
  constexpr auto tol = 0.0001220703;
  for(const auto& index:reconst_fcts){

    LogicalVector M_bool = is_na(m_Y(_,index));
    LogicalVector O_bool = !M_bool;
    GCV.set_bool(M_bool);
    
    if(alpha_nullable.isNull())
    {
      double max_bound = 0.0;
      for(size_t i = 0; i < r;++i){//r should be the nrow, ncol of m_cov -> diag is of length r
        max_bound += m_cov(i,i);
      }
      alpha_vec[column] = optimize(&GCV, std::numeric_limits<double>::epsilon(), max_bound * n, false, tol); //false -> minimization
    }else{
      alpha_vec[column] = alpha;
    }

    std::tuple<NumericVector, double, double, arma::vec, arma::uvec> resultKraus = reconstKraus_fun(m_Y, m_mean, m_cov, index, alpha_vec[column]);//forse dovrei cambiarla e passargli M_bool
  //reconstKraus["X_cent_reconst_vec"] is an arma::vec
    NumericVector X_reconst = std::get<0>(resultKraus);
    X_reconst_mat(_,column) = X_reconst + m_mean;
    df_vec[column] = std::get<1>(resultKraus);
    arma::uvec M_bool_arma = std::get<4>(resultKraus);
    arma::vec hi = std::get<3>(resultKraus);
    for(const auto& m:M_bool_arma)
    {
      W_reconst_mat(m, column) = 1 - hi[m - (r - hi.size())];
    }
    column++;
  }
  return(List::create(_["Y_reconst"] = X_reconst_mat,
                      _["alpha"]         = alpha_vec, 
                      _["df"]            = df_vec,
                      _["W_reconst"] = W_reconst_mat)); 

}


//extrapolation method for reconstruction from last observed period
List ReconstructionExtrapolation::reconstructCurve(Nullable<double> alpha = R_NilValue, bool all = FALSE, const Nullable<NumericVector>& t_points = R_NilValue, Nullable<int> K = R_NilValue, Nullable<int> maxBins = R_NilValue, Nullable<int> nRegGrid = R_NilValue) {
  int r;
  arma::vec periods;
  if(t_points.isNotNull()){
    periods = as<arma::vec>(t_points.get());
    r = periods.size();
  }else{
    stop("for extrapolation method you must provide a vector");
  }

  arma::mat m_Y_arma = as<arma::mat>(m_Y);
  int n = m_Y.ncol();
  arma::vec mean_slope(r-1); //mean row slope
  arma::uvec no_rec_fcts(n-m_reconst_fcts.size());
  //build vector of indeces not in m_reconst_fcts
  int row = 0;
  for( int i = 0; i < m_Y.ncol(); i++){
    bool found = std::binary_search(m_reconst_fcts.begin(),m_reconst_fcts.end(),i);//reconst_fcts is ordered
    if(!found){no_rec_fcts(row++) = i;}
  }

  const arma::mat& m_Y_no_rec = m_Y_arma.cols(no_rec_fcts);
  const arma::mat& last_row = m_Y_no_rec.row(r-1);
  const double& last_period = periods[r-1];
  for(int j = 0; j < r-1; j++)
  {
    arma::vec num = arma::conv_to<arma::vec>::from(last_row - m_Y_no_rec.row(j));
    mean_slope(j) = arma::mean(num)/(last_period - periods[j]);
  }

  arma::mat Y_reconstruct = m_Y_arma;
  for(const auto& index:m_reconst_fcts)
  {
    const arma::mat& col_Y_arma = m_Y_arma.col(index);
    const arma::uvec& id_no_na = arma::find_finite(col_Y_arma);
    const arma::uvec& id_na = arma::find_nonfinite(col_Y_arma);
    const arma::vec& observed_periods = periods(id_no_na);
    const double last_observed_period = observed_periods.back();
    const arma::vec& slopes = mean_slope(id_no_na);
    const double last_slope = slopes.back();
    const arma::vec& observed_col_Y = arma::conv_to<arma::vec>::from(col_Y_arma(id_no_na));
    const double last_observed_value = observed_col_Y.back();

    //reconstruct with equation of a straight line
    for(const auto& na:id_na)
    {
      Y_reconstruct(na,index) = last_observed_value + last_slope * (periods(na) - last_observed_period);
    }
  }
  return List::create(_["Y_reconst"] = wrap(Y_reconstruct),
                      _["mean_slopes"] = wrap(mean_slope));

}

std::pair<NumericMatrix,NumericVector> ReconstructionKL::smooth_cov(const NumericMatrix& Y_tilde,const NumericVector& Y_first,
                         int d, int i, int nbasis)
{
  NumericMatrix cov_mean(d, d);
  NumericMatrix cov_count(d,d);
  cov_count.fill(0.0);
  NumericMatrix cov_sum(d,d);
  cov_sum.fill(0.0);


  // Ccov.mean
  for (int idx = 0; idx < i; ++idx) {

    IntegerVector obs_points = seq_len(d) - 1;//0 based index
    obs_points = obs_points[!is_na((m_Y_preprocessed.second)(idx, _))];

    for(int i = 0; i < obs_points.size(); ++i) {
        for(int j = 0; j < obs_points.size(); ++j) {
            cov_count(obs_points[i], obs_points[j]) += 1;
        }
    }


    for(int i = 0; i < obs_points.size(); ++i)
    {
      for(int j = 0; j < obs_points.size(); ++j)
        cov_sum(obs_points[i],obs_points[j])  +=  Y_tilde(idx, obs_points[i]) * Y_tilde(idx, obs_points[j]);
    }
  }

  // G.0
  NumericMatrix G0(d, d);
  for (int row = 0; row < d; ++row) {
    for (int col = 0; col < d; ++col) {
      if (cov_count(row, col) == 0) {
          G0(row, col) = NA_REAL;
      }else {
          G0(row, col) = cov_sum(row, col) / cov_count(row, col);
      }
    }
  }

  NumericVector diag_G0(d);
  for (int i = 0; i < d; ++i) {
    diag_G0[i] = G0(i, i);
  }

  for (int i = 0; i < d; ++i) {
    G0(i, i) = NA_REAL;
  }


  // npc.0
  NumericVector row_vec(d * Y_first.size());
  int index = 0;
  for(const auto& y : Y_first) {
    for(int j = 0; j < d; ++j) {
      row_vec[index++] = y;
    }
  }

  NumericVector col_vec(Y_first.size() * d);
  for(int i = 0; i < d; ++i) {
    std::copy(Y_first.begin(), Y_first.end(), col_vec.begin() + i * Y_first.size());
  }

  /*Environment mgcv = Environment::namespace_env("mgcv");
  Function gam = mgcv["gam"];
  Function predict_gam = mgcv["predict.gam"];*/
  NumericVector g0(G0.begin(),G0.end());
  NumericVector weights(cov_count.begin(),cov_count.end());
  
  std::string formula = "G0 ~ te(row_vec, col_vec, k = " + std::to_string(nbasis) + ")";
  Formula f = Formula(formula);
  DataFrame data = DataFrame::create(_["G0"] = g0, _["row_vec"] = row_vec, _["col_vec"] = col_vec);
  List gamModel = gam(_["formula"] = f, _["data"] = data, _["weights"] = weights);//fit gam model
  DataFrame newdata = DataFrame::create(Named("row_vec") = row_vec, Named("col_vec") = col_vec);//data for prediction

  NumericVector predictions = predict_gam(gamModel, _["newdata"] = newdata);
  //reshape predictions vector into matrix
  NumericMatrix npc0(d, d);
  std::copy(predictions.begin(), predictions.end(), npc0.begin());
  NumericMatrix npc0_sym(d, d);

  // symmetrization
  for (int i = 0; i < d; ++i) {
    for (int j = i; j < d; ++j) {
          {
            npc0_sym(i, j) = (npc0(i, j) + npc0(j, i)) / 2.0;
            npc0_sym(j,i) = npc0_sym(i,j);
          }
    }
  }

  return std::make_pair(npc0_sym,diag_G0);
}

void ReconstructionKL::myfpca(std::vector<std::vector<double>>& Ly, const std::vector<std::vector<double>>& Lu, 
                                bool scores = false, bool center = true, Nullable<int> max_bins_nullable = R_NilValue, bool all = false){
  //int n = Ly.size();
  std::vector<size_t> id_vec = gen(Lu);
  std::vector<std::tuple<int, double, double>> ydata;
  ydata.reserve(id_vec.size());//giusto


  for(size_t i = 0; i < id_vec.size(); i++)
  {
    //id_vec 0-based per ora
    ydata.push_back(std::make_tuple(id_vec[i],unlist(Lu)[i],unlist(Ly)[i])); // .id,(time) .index, .value
  }

  constexpr int nbasis = 10;
  int max_bins = (max_bins_nullable.isNull())? 1000 : as<int>(max_bins_nullable); 
  //bool useSymm = false;
  //bool makePD = false;

  m_Y_preprocessed = irreg2mat(ydata, true, max_bins);//binning = true
  //argvals = m_Y_preprocessed.first

  IntegerVector reconst_fcts;
  if(all)//reconstruct all curves
  {
    reconst_fcts = seq_len(m_Y.ncol())-1;
  }else{
    reconst_fcts = (*this).reconst_fcts() - 1;
  }

  int d = m_Y_preprocessed.second.ncol();          //D 15
  int i = m_Y_preprocessed.second.nrow();          //I 5
  int i_pred = reconst_fcts.size();                //I.PRED

  NumericMatrix Y_pred(i_pred, d);  //m_Y_preprocessed.second è la matrice
  int row = 0;
  for(const auto& index: reconst_fcts)
  {
    Y_pred(row++,_) = m_Y_preprocessed.second(index,_);
  }

  std::vector<std::vector<double>> observed_period; //argvalsO
  observed_period.reserve(i_pred);

  for(int ii = 0; ii < i_pred; ++ii) //reconst_fcts.size() = Y_pred.nrow()
  {
    std::vector<double> obs;
    obs.reserve(d);
    for(int j = 0; j < d; ++j)
    {
      if(!NumericVector::is_na(Y_pred(ii,j)))
        obs.push_back(m_Y_preprocessed.first[j]);
    }
    obs.shrink_to_fit();
    observed_period.push_back(obs);
  }

  m_observed_period = observed_period;
      
  std::vector<double> d_vec;
  d_vec.reserve(i * d); 
  for (const double& val : m_Y_preprocessed.first) {
    for (int j = 0; j < i; ++j) {//m_Y_preprocessed.first = argvals = t.points  
        d_vec.push_back(val);        
    }
  }
  d_vec.shrink_to_fit();//non servirebbe

  std::vector<int> id;
  id.reserve(i * d); 
  for (int j = 0; j < i; ++j) {
      for (int dd = 0; dd < d; ++dd) {
          id.push_back(j);
      }
  }
  id.shrink_to_fit();

  NumericMatrix Y_tilde(clone(m_Y_preprocessed.second)); //cercare se da altre parti avrei dovuto usare clone per la shallow copy
  NumericVector yfirst((m_Y_preprocessed.first).begin(), (m_Y_preprocessed.first).end());   //t_points
  NumericVector mu;
  if(center)
  {
    NumericMatrix Y_mat = m_Y_preprocessed.second;
    NumericVector vec(Y_mat.begin(), Y_mat.end());

    NumericVector dvec(d_vec.begin(), d_vec.end());

    DataFrame data = DataFrame::create(_["Y"] = vec, _["d.vec"] = dvec);
    std::string formula_str = "Y ~ s(d.vec, k = " + std::to_string(nbasis) + ")";
    Formula f = Formula(formula_str);

    List gam0 = gam(_["formula"] = f, _["data"] = data);

    DataFrame newdata = DataFrame::create(_["d.vec"] = yfirst);
    mu = predict_gam(gam0, _["newdata"] = newdata);

    
    for (int row = 0; row < i; ++row) {
      for(int col = 0; col < d; ++col){
        Y_tilde(row,col) -= mu[col];
      }
    }
  }else{
    mu = NumericVector(d,0.0);
  }
  m_mu = mu;//giusto
  //problema principale -> t_points è Y.first. potrei non usare Y.first
  //cov
  std::pair<NumericMatrix, NumericVector> cov_smooth = smooth_cov(Y_tilde, yfirst, d, i, nbasis);//if(!useSymm), perchè setta useSymm = FALSE
  NumericMatrix npc0 = cov_smooth.first; //giusto
  NumericVector diagG0 = cov_smooth.second; //giusto

  //numerical integration for calculation of eigenvalues (see Ramsay & Silverman, Ch. 8)
  std::tuple<std::vector<arma::vec>, std::vector<double>, List, arma::mat, std::vector<arma::mat>, std::vector<NumericMatrix>, 
           arma::vec,std::vector<arma::vec>, std::vector<arma::vec>, std::vector<arma::uvec>, double, arma::mat>  ret = eigen(m_Y_preprocessed.first, observed_period, npc0, pev, Y_pred, mu, diagG0, true, reconst_fcts);
  
  //fissa i nuovi data member
  m_muO = std::get<0>(ret);//giusto
  m_scoresO = std::get<1>(ret);//no?
  m_CE_scoresO = std::get<2>(ret);//giusto
  m_efunctions = std::get<3>(ret);//giusto
  m_efunctionsO = std::get<4>(ret);//giusto
  m_efun_reconst = std::get<5>(ret);//giusto
  m_evalues = std::get<6>(ret);//no?
  m_evaluesOO = std::get<7>(ret);//giusto
  m_obs_argvalsO = std::get<8>(ret);//giusto
  m_locO = std::get<9>(ret);//controlla se ha senso tenerlo e farlo restituire, o se non ho già questa informazione da altre parti  
  m_sigma2 = std::get<10>(ret);//giusto
  m_cov_est = std::get<11>(ret);//giusto
}


List ReconstructionKL::reconstructCurve(Nullable<double> alpha = R_NilValue, bool all = FALSE, const Nullable<NumericVector>& t_points = R_NilValue, Nullable<int> K = R_NilValue, Nullable<int> maxBins = R_NilValue, Nullable<int> nRegGrid_nullable = R_NilValue)
{ 
  NumericVector t_points_;
  if(t_points.isNull()){
    stop("you must provide a vector of t.points");
  }else{
    t_points_ = t_points.get();
  }
  int n = m_Y.ncol();
  int r = m_Y.nrow();
  std::vector<std::vector<double>> Y_list(n);
  std::vector<std::vector<double>> U_list(n);
  NumericVector reconst_fcts;
  if(all)
  {
    reconst_fcts = seq_len(n) - 1;
  }else{
    reconst_fcts = (*this).reconst_fcts() - 1;
  }

  for(int i = 0; i < n; i++)
  {
    std::vector<double> keep;
    std::vector<double> keep_t;
    keep.reserve(r);
    keep_t.reserve(r);
    if(std::binary_search(m_reconst_fcts.begin(),m_reconst_fcts.end(), i))
    {
      for(int j = 0; j < r; j++)
      {
        if(!NumericVector::is_na(m_Y(j,i)))
          {
            keep.push_back(m_Y(j,i));
            keep_t.push_back(t_points_[j]);
          }
      }
      keep.shrink_to_fit();
      keep_t.shrink_to_fit();

    }else{
      NumericVector col = m_Y(_,i);
      keep = as<std::vector<double>>(col);
      keep_t = as<std::vector<double>>(t_points_);
    }
    U_list.push_back(keep_t);
    Y_list.push_back(keep);
  }//corretto
  //created Y_list and U_list 
  //myfpca setta dei data member. reonstructCurve() non può essere const
  myfpca(Y_list, U_list, false, true, maxBins, all);//all = false vuol dire che non deve ricostruire tutto

  m_length_reconst_fcts = reconst_fcts.size();
  std::vector<double> K_vec;
  K_vec.reserve(m_length_reconst_fcts);
  std::vector<arma::vec> Y_reconstr_list(m_length_reconst_fcts);
  std::vector<arma::vec> W_reconst_list(m_length_reconst_fcts);
  std::vector<arma::vec> U_reconst_list(m_length_reconst_fcts);
  NumericVector x((m_Y_preprocessed.first).begin(), (m_Y_preprocessed.first).end());
  for(size_t i = 0; i < m_length_reconst_fcts; ++i){
    const arma::vec& obs = m_obs_argvalsO[i];
    const std::vector<double>& argvalsO_i = m_observed_period[i];//argvalsO[i] in R
    if(!(obs[0] == argvalsO_i[0] && obs[obs.size()-1] == argvalsO_i[argvalsO_i.size()-1]))
    {
      stop("The range of obs_argvalsO of the fragment must equal the range of argvalsO");
    }
    
    if(K.isNull())
    {
      K_vec.push_back(gcvKneipLiebl(m_mu, m_Y_preprocessed, argvalsO_i, m_locO[i], m_cov_est, m_sigma2, m_method, pev));//giusto
    }else{K_vec.push_back(as<int>(K));}

    NumericVector fragmO_presmooth = NumericVector();
    if(m_method == "KLAl"){
        NumericVector y = (m_Y_preprocessed.second)(reconst_fcts[i],_);
        LogicalVector no_na = !is_na(y);
        NumericVector y_c = y[no_na];
        NumericVector obs_argvals(obs.begin(), obs.end());
        NumericVector argvalsO_i_vector(argvalsO_i.begin(), argvalsO_i.end());

        List smooth_fit = smooth_spline(_["y"] = y_c, _["x"] = obs_argvals);
        List fragmO_presmooth_list = predict(smooth_fit, _["x"] = argvalsO_i_vector);
        fragmO_presmooth = fragmO_presmooth_list["y"];
    }//giusto

    //fragmO_presmooth empty if method == "KLNoAl"
    std::pair<arma::vec, arma::vec> result = reconstKL_fun(m_mu, m_Y_preprocessed.first, m_locO[i], m_CE_scoresO[i], 
                                                            m_efun_reconst[i], fragmO_presmooth, K_vec[i],
                                                            m_evaluesOO[i], m_observed_period[i], m_cov_est);
    Y_reconstr_list[i] = result.first;
    W_reconst_list[i] = result.second;
    U_reconst_list[i] = x;
  }
  if(nRegGrid_nullable.isNotNull())
  {
    int nRegGrid = as<int>(nRegGrid_nullable);
    double start = m_Y_preprocessed.first[0];
    double end = m_Y_preprocessed.first[m_Y_preprocessed.first.size()-1];
    double step = (end - start)/(nRegGrid - 1);
    NumericVector sequence(nRegGrid);
    for (int i = 0; i < nRegGrid; ++i) {
      sequence[i] = start + i * step;
    }
    sequence = sequence * step + start;
  
    for(size_t i = 0; i < m_length_reconst_fcts; ++i)
    {
      NumericVector y((Y_reconstr_list[i]).begin(), (Y_reconstr_list[i]).end());
      List reconstr_on_grid = spline(_["y"] = y, _["x"] = x, _["xout"] = sequence);
      Y_reconstr_list[i] = as<arma::vec>(reconstr_on_grid["y"]);
      U_reconst_list[i] = as<arma::vec>(reconstr_on_grid["x"]);
    }
  }
  NumericVector K_vec_(K_vec.begin(),K_vec.end());
  return List::create(_["Y_reconst_list"] = Y_reconstr_list, _["U_reconst_list"] = U_reconst_list, 
                      _["W_reconst_list"] = W_reconst_list, _["K"] = K_vec_);
}
