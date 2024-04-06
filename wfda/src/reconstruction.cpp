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


List ReconstructionKraus::reconstructCurve(Nullable<double> alpha_nullable = R_NilValue, bool all = FALSE, const Nullable<NumericVector>& t_points = R_NilValue, Nullable<int> K = R_NilValue, Nullable<int> maxBins = R_NilValue, Nullable<int> nRegGrid = R_NilValue) {
//dovrei avere già mean_vec e cov_mat nella classe appena chiamo il costruttore
  double alpha;
  if(alpha_nullable.isNotNull()){alpha = as<double>(alpha_nullable);} 
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
  std::fill(W_reconst_mat.begin(),W_reconst_mat.end(),1);//posso fare così grazie a come sono salvate le NumericMatrix in memoria

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
    
    if(alpha_nullable.isNull())//R_NilValue = NULL in R, in attesa di capire come gestire NULL, metto 0.0 di default
    {
      double max_bound = 0.0;
      for(int i = 0; i < r;++i){//r should be the nrow, ncol of m_cov -> diag is of length r
        max_bound += m_cov(i,i);
      }
      alpha_vec[column] = optimize(&GCV, std::numeric_limits<double>::epsilon(), max_bound, false, 1e-3); //false -> minimization
    }else{
      alpha_vec[column] = alpha;
    }

    std::tuple<NumericVector, double, double, arma::vec, arma::uvec> resultKraus = reconstKraus_fun(m_Y, m_mean, m_cov, index, alpha_vec[column]);//forse dovrei cambiarla e passargli M_bool
  //reconstKraus["X_cent_reconst_vec"] is an arma::vec
    NumericVector X_reconst = std::get<0>(resultKraus);
    X_reconst_mat(_,column) = X_reconst + m_mean;
    df_vec[column] = std::get<1>(resultKraus);//check
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
  arma::vec mean_slope(r-1); //slope media per ogni riga
  arma::uvec no_rec_fcts(n-m_reconst_fcts.size());
  //build vector of indeces not in m_reconst_fcts
  int row = 0;
  for( int i = 0; i < m_Y.ncol();i++){
    bool found = std::binary_search(m_reconst_fcts.begin(),m_reconst_fcts.end(),i);//works cause reconst_fcts is an ordered vector
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

  constexpr int nbasis = 10;//vedi se ha senso mettere qua constexpr
  int max_bins = (max_bins_nullable.isNull())? 1000 : as<int>(max_bins_nullable); 
  //bool useSymm = false;//se sono sempre false finsico sempre negli stessi branch
  //bool makePD = false;

  std::pair<std::vector<double>,NumericMatrix> Y = irreg2mat(ydata, true, max_bins);//binning = true
  m_Y_preprocessed = Y;
  //argvals è Y.first
  IntegerVector reconst_fcts;
  if(all)//reconstruct all curves
  {
    reconst_fcts = seq_len(m_Y.ncol())-1;//lo devo correggere anche in Kraus e mettere -1
  }else{
    reconst_fcts = wrap(m_reconst_fcts);
  }

  int d = Y.second.ncol();          //D 15
  int i = Y.second.nrow();          //I 5
  int i_pred = reconst_fcts.size(); //I.PRED

  NumericMatrix Y_pred(i_pred, d);  //Y.second è la matrice
  int row = 0;
  for(const auto& index: reconst_fcts)
  {
    Y_pred(row++,_) = Y.second(index,_);
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
        obs.push_back(Y.first[j]);
    }
    obs.shrink_to_fit();
    observed_period.push_back(obs);
  }

  m_observed_period = observed_period;

  //argvals = Y.first
      
  std::vector<double> d_vec;
  d_vec.reserve(i * d); 
  for (const double& val : Y.first) {
    for (int j = 0; j < i; ++j) {//Y.first = argvals = t.points .....  
        d_vec.push_back(val);        
    }
  }
  d_vec.shrink_to_fit();//non servirebbe

  std::vector<int> id;
  id.reserve(i * d); 
  for (int j = 0; j < i; ++j) {//controlla che non sia da j = 1 a j = i, ma solo perchè devo capire se servono gli indici di R o C++, penso C++
      for (int dd = 0; dd < d; ++dd) {
          id.push_back(j);
      }
  }
  id.shrink_to_fit();

  NumericMatrix Y_tilde(clone(Y.second)); //cercare se da altre parti avrei dovuto usare clone per la shallow copy
  NumericVector yfirst = wrap(Y.first);//t_points
  NumericVector mu;
  if(center)
  {
    NumericMatrix Y_mat = Y.second;
    NumericVector vec(Y_mat.begin(), Y_mat.end());
    Environment mgcv = Environment::namespace_env("mgcv");    
    Function gam = mgcv["gam"];
    Function predict_gam = mgcv["predict.gam"];

    NumericVector dvec = wrap(d_vec); //forse posso passare direttamente d_vec e chiamato in automatico wrap

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
  m_mu = mu;
  //problema principale -> t_points è Y.first. potrei non usare Y.first
  //cov
  std::pair<NumericMatrix, NumericVector> cov_smooth = smooth_cov(Y.second, Y_tilde, yfirst, d, i, nbasis);//if(!useSymm), perchè setta useSymm = FALSE
  NumericMatrix npc0 = cov_smooth.first; //giusto
  NumericVector diagG0 = cov_smooth.second; //giusto

  //numerical integration for calculation of eigenvalues (see Ramsay & Silverman, Ch. 8)
  std::tuple<std::vector<arma::vec>, std::vector<double>, List, arma::mat, std::vector<arma::mat>, std::vector<NumericMatrix>, 
           arma::vec,std::vector<arma::vec>, std::vector<arma::vec>, std::vector<arma::uvec>, double, arma::mat>  ret = eigen(Y.first, observed_period, npc0, 0.99, Y_pred, mu, diagG0, true, reconst_fcts);
  
  //fissa i nuovi data member
  
  m_muO = std::get<0>(ret);
  m_scoresO = std::get<1>(ret);
  m_CE_scoresO = std::get<2>(ret);
  m_efunctions = std::get<3>(ret);
  m_efunctionsO = std::get<4>(ret);
  m_efun_reconst = std::get<5>(ret);
  m_evalues = std::get<6>(ret);
  m_evaluesOO = std::get<7>(ret);
  m_obs_argvalsO = std::get<8>(ret);
  m_locO = std::get<9>(ret);
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
  if(all){reconst_fcts = seq_len(n) - 1;}else{reconst_fcts = wrap(m_reconst_fcts);}
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

  size_t length_reconst_fcts = m_obs_argvalsO.size();//per come è stata costruita in eigen è di dimensione pari a lunghezza reconst_fcts
  std::vector<double> K_vec;
  K_vec.reserve(length_reconst_fcts);
  std::vector<arma::vec> Y_reconstr_list(length_reconst_fcts);
  std::vector<arma::vec> W_reconst_list(length_reconst_fcts);
  std::vector<arma::vec> U_reconst_list(length_reconst_fcts);
  Environment stats = Environment::namespace_env("stats");
  NumericVector x = wrap(m_Y_preprocessed.first);
  for(size_t i = 0; i < length_reconst_fcts; ++i){
    arma::vec obs = m_obs_argvalsO[i];
    std::vector<double> argvalsO_i = m_observed_period[i];//argvalsO[i] in R
    if(!(obs[0] == argvalsO_i[0] && obs[obs.size()-1] == argvalsO_i[argvalsO_i.size()-1]))
    {
      stop("The range of obs_argvalsO of the fragment must equal the range of argvalsO");
    }
    
    if(K.isNull())//per ora ho dato il default a 0, ma poi mettere Nullable
    {
      //std::string method = "KLAl4";//pev = 0.99
      K_vec.push_back(gcvKneipLiebl(m_mu, m_Y_preprocessed, argvalsO_i, m_locO[i], m_cov_est, m_sigma2, m_method, 0.99));
    }else{K_vec.push_back(as<int>(K));}

    NumericVector fragmO_presmooth = NumericVector();
    if(m_method == "KLAl"){
        NumericVector y = (m_Y_preprocessed.second)(reconst_fcts[i],_);
        LogicalVector no_na = !is_na(y);
        NumericVector y_c = y[no_na];
        NumericVector obs_argvals = wrap(m_obs_argvalsO[i]); //wrap the arma::vec
        NumericVector argvalsO_i_vector = wrap(argvalsO_i);

        Function smooth_spline = stats["smooth.spline"];
        List smooth_fit = smooth_spline(_["y"] = y_c, _["x"] = obs_argvals);

        Function predict = stats["predict"];
        List fragmO_presmooth_list = predict(smooth_fit, _["x"] = argvalsO_i_vector);
        fragmO_presmooth = fragmO_presmooth_list["y"];
    }

    //fragmO_presmooth empty if method == "KLNoAl"
    std::tuple<arma::vec, arma::vec> result = reconstKL_fun(m_mu, m_Y_preprocessed.first, m_locO[i], m_CE_scoresO[i], 
                                                            m_efun_reconst[i], fragmO_presmooth, K_vec[i],
                                                            m_evaluesOO[i], m_observed_period[i], m_cov_est);
    Y_reconstr_list[i] = std::get<0>(result);
    W_reconst_list[i] =std::get<1>(result);
    U_reconst_list[i] = x;
  }
  if(nRegGrid_nullable.isNotNull())
  {
    int nRegGrid = as<int>(nRegGrid_nullable);
    Function spline = stats["spline"];
    double start = m_Y_preprocessed.first[0];
    double end = m_Y_preprocessed.first[m_Y_preprocessed.first.size()-1];
    double step = (end - start)/(nRegGrid - 1);
    NumericVector sequence(nRegGrid);
    for (int i = 0; i < nRegGrid; ++i) {
      sequence[i] = start + i * step;
    }
    sequence = sequence * step + start;
  
    for(size_t i = 0; i < length_reconst_fcts; ++i)
    {
      NumericVector y = wrap(Y_reconstr_list[i]);
      List reconstr_on_grid = spline(_["y"] = y, _["x"] = x, _["xout"] = sequence);
      Y_reconstr_list[i] = as<arma::vec>(reconstr_on_grid["y"]);
      U_reconst_list[i] = as<arma::vec>(reconstr_on_grid["x"]);
    }
  }
  NumericVector K_vec_ = wrap(K_vec);
  return List::create(_["Y_reconst_list"] = Y_reconstr_list, _["U_reconst_list"] = U_reconst_list, 
                      _["W_reconst_list"] = W_reconst_list, _["K"] = K_vec_);
}