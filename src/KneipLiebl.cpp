#include <RcppArmadillo.h>
#include "helper_functions.h"
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
std::pair<arma::vec, arma::vec> 
reconstKL_fun(const NumericVector& mu, const std::vector<double>& argvals, const arma::uvec& locO, 
              const arma::vec& scoresO, const NumericMatrix& efunc_r, const NumericVector& fragmO, int k,
              const arma::vec& evaluesO = arma::vec(), const std::vector<double>& argvalsO = std::vector<double>(), const arma::mat& cov = arma::mat())
{
  int K = std::min(k, efunc_r.ncol());
  arma::mat efunc_r_arma = as<arma::mat>(efunc_r);
  arma::vec mu_arma = as<arma::vec>(mu);//giusto
  arma::mat extracted_scores = arma::reshape(scoresO.subvec(0, K - 1), 1, K);
  arma::vec reconstr = mu_arma + ((efunc_r_arma.cols(0,K-1)) * extracted_scores.t());//giusto
  //                   15 x 1  +  15 x k                     *          k x 1 = 15 x 1 + 15 x 1

  if(fragmO.size() != 0)//metodo4
  {
    //arma::uword min = locO.min();//locO min sarà sempre 0?
    //reconstr.subvec(0, min) += fragmO[0] - reconstr[min];
    arma::uword max = locO.max();
    reconstr.subvec(max + 1, argvals.size() - 1) += fragmO[fragmO.size() - 1] - reconstr[max];
    reconstr(locO) = as<arma::vec>(fragmO);//giusto
  }
  arma::vec weights_reconst;
  if(cov.is_empty() || evaluesO.is_empty())
  {
    weights_reconst = arma::vec();
  }else{
    arma::uvec locM = arma::regspace<arma::uvec>(locO.max() + 1, argvals.size() - 1);
    arma::vec diag_cov = arma::diagvec(cov);
    arma::mat efunc_r_submat = efunc_r_arma.cols(0, K-1) % efunc_r_arma.cols(0, K-1);//element wise
    arma::vec e = efunc_r_submat * evaluesO.subvec(0, K-1) ;
    arma::vec v2_reconstr = diag_cov.elem(locM) - e.elem(locM);
    arma::vec v_hat_reconstr(v2_reconstr.size());
    if(arma::all(v2_reconstr > 1e-14)){
      arma::vec v_reconstr = arma::sqrt(v2_reconstr);
      double h0 = 0.2*arma::max(v_reconstr);
      std::transform(v_reconstr.begin(), v_reconstr.end(), v_hat_reconstr.begin(),
                     [h0](double value){ return (h0 < value ? value :h0);});

    }else{
      v_hat_reconstr.fill(0.01);
    }
    weights_reconst.resize(argvals.size());
    weights_reconst.fill(1);
    weights_reconst.elem(locM) -= v_hat_reconstr / arma::sqrt(diag_cov.elem(locM));
  }
  return std::make_pair(reconstr, weights_reconst);
}


std::vector<std::tuple<int, double, double>> //a me viene che ydata è già complete
find_complete_tuple(const std::vector<std::tuple<int, double, double>>& y_data)
{
  std::vector<std::tuple<int, double, double>> y_data_complete;
  y_data_complete.reserve(y_data.size()); //maximum

  for(const auto& t:y_data)//vedere dove altro mettere const nei range for (penso dappertutto)
  {
    if(!NumericVector::is_na(std::get<1>(t)) && !NumericVector::is_na(std::get<2>(t)))//forse basta controllare uno dei due
    {
      y_data_complete.push_back(t);
    }
  }
  y_data_complete.shrink_to_fit();
  return y_data_complete;
}


//irreg2mat <- function(ydata, binning = FALSE, maxbins = 1000)
std::pair<std::vector<double>,NumericMatrix> irreg2mat(const std::vector<std::tuple<int, double, double>>& y_data, 
                                                       bool binning, int max_bins){

  //crea una copia
  std::vector<std::tuple<int, double, double>> y_data_complete(find_complete_tuple(y_data));
  std::set<double> bins; 

  std::vector<int> ids;
  ids.reserve(y_data_complete.size());
  for(const auto& t: y_data_complete)
  {
    ids.push_back(std::get<0>(t));//newid. deve essere il vettore di indici di riga di Y
    bins.insert(std::get<1>(t));
  }//bins è per forza uguale a t_points

  //ids already ordered, no sorting beforehand needed
  std::vector<int> ids_copy = ids;//copy since std::unique modifies vector
  int nobs = std::distance(ids_copy.begin(),std::unique(ids_copy.begin(),ids_copy.end()));//corretto
  bool condition = binning && bins.size() > static_cast<size_t>(max_bins);//false
  std::vector<double> binvalues = make_bins(bins, max_bins, condition);//bins sarà modificato, chiamo make_bins con la reference

  //associa a std::get<1>(t) un bin => crea classi di .index
  //estrai vettore dal secondo elemento della tupla
  std::vector<double> index_vector;
  index_vector.reserve(y_data_complete.size());
  for(const auto& t:y_data_complete)
  {
    index_vector.push_back(std::get<1>(t));
  }
  index_vector.shrink_to_fit();

  std::vector<std::pair<double, double>> classes = cut(index_vector, bins);//this is newindex in R
  //associate each class to a column index
  std::set<std::pair<double,double>> sorted_classes(classes.begin(),classes.end());
  //associate index
  std::vector<int> column_indices(classes.size()); //vedi se trasformare in size_t       

  for (size_t i = 0; i < classes.size(); ++i) {
      auto it = std::find(sorted_classes.begin(), sorted_classes.end(), classes[i]);//lo trova per forza. find su un vettore restituisce iteratore
      column_indices[i] = std::distance(sorted_classes.begin(), it);
  }
 
  //inizializza ad NA_REAL, alcune curve non hanno ossservazioni in dati bins
  NumericMatrix Y(nobs, bins.size()-1);
  Y.fill(NA_REAL);
  //nobs*(bins.size()-1) >= y_data_complete.size()
  for(size_t i = 0; i < y_data_complete.size();++i){
    Y(ids[i],column_indices[i]) = (std::get<2>(y_data_complete[i]));//ydata$.value
  }

  std::pair<std::vector<double>,NumericMatrix> Y_bins = std::make_pair(binvalues,Y);//binvalues is colnames(Y)
  return Y_bins;
  
}


std::vector<double> quadWeights(const std::vector<double>& argvals, const std::string& method = "trapezoidal") {
    int D = argvals.size();
    std::vector<double> ret(D, 0.0);

    if (method == "trapezoidal") {
        ret[0] = 0.5 * (argvals[1] - argvals[0]);
        for (int i = 1; i < D - 1; i++) {
            ret[i] = 0.5 * (argvals[i + 1] - argvals[i - 1]);
        }
        ret[D - 1] = 0.5 * (argvals[D - 1] - argvals[D - 2]);
    } else if (method == "midpoint") {
        ret[0] = 0.0; 
        for (int i = 1; i < D; i++) {
            ret[i] = argvals[i] - argvals[i - 1];
        }
    } else {
        stop("function quadWeights: choose either trapezoidal or midpoint quadrature rule");
    }

    return ret;
}


double weighted_mean(const std::vector<double>& v, const std::vector<double>& w)
{
  if(v.size() != w.size())
    stop("vector of weights must be the same size as of input vector");
  double sum = 0.0;
  double sum_w = 0.0;
  // #pragma omp parallel for reduction(+:sum, sum_w)
  for(size_t i = 0; i < v.size(); ++i)
  {
    sum += v[i]*w[i];
    sum_w += w[i];//this does not necessarily sum to 1 ?
  }
  return sum/sum_w; 
}


double trapezioidal_rule(const arma::vec& values, const arma::vec& knots){
  double integral = 0.0;
  size_t n = knots.size();

  for (size_t i = 1; i < n; ++i) {
      integral += (knots[i] - knots[i - 1]) * (values[i] + values[i - 1]) / 2.0;
  }

  return integral;
}

std::tuple<std::vector<arma::vec>, std::vector<double>, List, arma::mat, std::vector<arma::mat>, std::vector<NumericMatrix>, 
           arma::vec,std::vector<arma::vec>, std::vector<arma::vec>, std::vector<arma::uvec>, double, arma::mat> 
           eigen(const std::vector<double>& argvals, const std::vector<std::vector<double>>& argvalsO,
           const NumericMatrix& npc_0, 
           double pev, const NumericMatrix& Y_pred, 
           const NumericVector& mu, const NumericVector& diagG0, bool CEScores, 
           const IntegerVector& reconst_fcts) {
  
  // numerical integration for calculation of eigenvalues
  std::vector<double> w = quadWeights(argvals);//default is trapezioidal
  arma::vec w_arma = arma::conv_to<arma::vec>::from(w);
  arma::mat Wsqrt = arma::diagmat(arma::sqrt(w_arma));
  arma::mat Winvsqrt = arma::diagmat(1 / arma::sqrt(w_arma));
  arma::mat V = Wsqrt * as<arma::mat>(npc_0) * Wsqrt;

  arma::vec evalues;
  arma::mat eigenvectors;
  arma::eig_sym(evalues, eigenvectors, V);//ascending orderered eigenvalues
  arma::uvec sort_index = arma::sort_index(evalues, "desc");
  //descending order, coherence with R
  evalues = evalues(sort_index);
  eigenvectors = eigenvectors.cols(sort_index);
  evalues.transform([](double val) { return val <= 0 ? 0.0 : val; });//eigenvalues >= 0
  int npc = arma::sum(evalues > 0) - 1;//position index: subtract -1
  double sum_pos = arma::sum(evalues(arma::find(evalues > 0)));

  if (!NumericVector::is_na(pev)) { 
    arma::uword npc_index;

    if (sum_pos > 0) {
      arma::vec cumsum_pos = arma::cumsum(evalues(arma::find(evalues > 0))); 
      double threshold = pev * sum_pos;
      arma::uvec indices = arma::find(cumsum_pos >= threshold);

      if(!indices.is_empty())
        {npc_index = indices[0]; npc = static_cast<int>(npc_index);}
    }
  }

  //only select eigenvectors and corresponding eigenvalues up to npc-th
  arma::mat selected_eigenvectors = eigenvectors.cols(0, npc);
  arma::mat efunctions = Winvsqrt * selected_eigenvectors;
  arma::vec selected_evalues = evalues.subvec(0, npc);

  
  //estimated covariance function
  arma::mat cov_est = efunctions * arma::diagmat(selected_evalues) * efunctions.t();
  //numerical integration for estimation of sigma2
  double T_len = argvals[argvals.size() - 1] - argvals[0];//interval length

  double lower_bound = argvals[0] + 0.25 * T_len;
  auto it = std::lower_bound(argvals.begin(), argvals.end(), lower_bound);//use it on ordered container
  int T1_min = std::distance(argvals.begin(), it);

  double upper_bound = argvals[argvals.size() - 1] - 0.25 * T_len;
  auto it_max = std::upper_bound(argvals.begin(), argvals.end(), upper_bound);
  if(it_max != argvals.begin())
  {
    --it_max;
  }

  int T1_max = std::distance(argvals.begin(), it_max);// = 0 if all > upper_bound
  arma::vec diag_diff = as<arma::vec>(diagG0) - arma::diagvec(cov_est);
  arma::vec sub_diag = diag_diff.subvec(T1_min, T1_max);
  std::vector<double> argvals_subset(argvals.begin() + T1_min, argvals.begin() + T1_max + 1);

  std::vector<double> w2 = quadWeights(argvals_subset);//default: trapezioidal

  //sigma
  std::vector<double> sub_diag_vec(sub_diag.begin(), sub_diag.end());
  double sigma2 = std::max(weighted_mean(sub_diag_vec, w2), 0.0);
  //computations for observed fragments
  size_t length_reconst_fcts = reconst_fcts.size();
  std::vector<arma::vec> muO(length_reconst_fcts), obs_argvalsO(length_reconst_fcts), evaluesOO(length_reconst_fcts);
  std::vector<double> scoresO(length_reconst_fcts);
  std::vector<arma::mat> efunctionsO(length_reconst_fcts);
  std::vector<NumericMatrix> efun_reconst(length_reconst_fcts);
  std::vector<arma::uvec> locOO(length_reconst_fcts);
  List CE_scoresO(length_reconst_fcts);
  
  
  for (size_t i = 0; i < length_reconst_fcts; ++i) {
    // Numerical integration for calculation of eigenvalues
    std::vector<double> w = quadWeights(argvalsO[i]);
    arma::vec w_arma = arma::conv_to<arma::vec>::from(w);
    arma::mat Wsqrt = arma::diagmat(arma::sqrt(w_arma));
    arma::mat Winvsqrt = arma::diagmat(1 / arma::sqrt(w_arma));

    arma::uvec locO(argvalsO[i].size());//0-indexed
    for(arma::uword j = 0; j < locO.size(); ++j)
    {
      locO[j] = j; //match() su R
    }
    locOO[i] = locO;

    // CovOO
    arma::mat VO = Wsqrt * cov_est.submat(locO,locO) * Wsqrt;
    arma::vec evaluesO;
    arma::mat eigenvectorsO;
    arma::eig_sym(evaluesO, eigenvectorsO, VO);
    arma::uvec sort_indices = arma::sort_index(evaluesO, "desc");
    evaluesO = evaluesO(sort_indices);
    eigenvectorsO = eigenvectorsO.cols(sort_indices);
    evaluesO.transform([](double val){ return val <= 0 ? 0.0 : val;});
    size_t npcO = arma::sum(evaluesO > 0);
    if(npcO == 0){stop("no eigen values greater than 0");}else{npcO--;}

    double sum_pos_O = arma::sum(evaluesO(arma::find(evaluesO > 0)));

    if (!NumericVector::is_na(pev)) { //riguardare default di pev, in R è NULL, qua sto facendo come se fosse NA_REAL
      arma::uword npc_index_i;

      if (sum_pos > 0) {
        arma::vec cumsum_pos_O = arma::cumsum(evaluesO(arma::find(evaluesO > 0)));
        double threshold = pev * sum_pos_O;
        arma::uvec indices = arma::find(cumsum_pos_O >= threshold);

        if(!indices.is_empty())//in realtà non dovrebbe essere mai empty
          {npc_index_i = indices[0]; npcO = static_cast<int>(npc_index_i);}
      }
    }

    arma::mat efunctionsO_i = Winvsqrt * eigenvectorsO.cols(0, npcO);
    arma::vec evaluesO_i = evaluesO.subvec(0, npcO); //giusto

    arma::mat D_inv = arma::diagmat(1/evaluesO_i);
    arma::mat Z = efunctionsO_i;
    NumericVector Y_cent = Y_pred(i,_) - mu;
    arma::uvec obs_locO(Y_pred.ncol() - sum(is_na(Y_cent)));
    for(arma::uword j = 0; j < obs_locO.size(); ++j){
      obs_locO[j] = j;
    }
    arma::vec obs_argvalsO_i = (arma::conv_to<arma::vec>::from(argvalsO[i])).elem(obs_locO);
    obs_argvalsO[i] = obs_argvalsO_i;
    arma::vec Y_cent_arma = as<arma::vec>(Y_cent);
    arma::uvec no_na = arma::find_finite(Y_cent_arma);
    // Calculate CE_scoresO
    bool size = false;
    if (CEScores) {
      if (sigma2 == 0) {
        sigma2 = 1e-6;
      }
      if (obs_locO.size() < npcO + 1) {//+1 cause it was an index and now a size
        npcO = obs_locO.size();
        size = true;//if here it is a size not an index
      }
      arma::mat Zcur = Z.submat(arma::span(obs_locO[0],obs_locO[obs_locO.size()-1]),arma::span(0, npcO - 1*CEScores*size));
      arma::mat ZtZ_sD_inv = arma::inv(Zcur.t() * Zcur + sigma2 * D_inv.submat(arma::span(0, npcO - 1*CEScores*size),arma::span(0, npcO - 1*CEScores*size)));
      arma::vec CE_scoresO_i = ZtZ_sD_inv * Zcur.t() * Y_cent_arma(no_na);//arma::vec column vectors
      CE_scoresO[i] = CE_scoresO_i;
    } else {
      CE_scoresO[i] = NA_REAL;
    }

    arma::mat efunctionsO_i_sub = efunctionsO_i.rows(obs_locO);
    for(arma::uword j = 0; j < efunctionsO_i_sub.n_cols;++j)
    {
      arma::vec column = efunctionsO_i_sub.col(j);
      double integral = trapezioidal_rule(column % Y_cent_arma(no_na), obs_argvalsO_i);
      scoresO[i] = integral;
    }

    //reconstructive eigenfunctions
    NumericMatrix efun_reconst_i(argvals.size(), npcO + 1*(!size));//giusto
    efun_reconst_i.fill(NA_REAL);     
    arma::mat rows = cov_est.rows(locO);
    for(int k = 0; k < efun_reconst_i.ncol(); ++k)
    {
      for(int r = 0; r < efun_reconst_i.nrow(); ++r){
          arma::vec c = rows.col(r);
          double integral = trapezioidal_rule(efunctionsO_i.col(k) % c, arma::conv_to<arma::vec>::from(argvalsO[i]));
          efun_reconst_i(r,k) = integral/evaluesO_i[k];
      }
    }
    
    efun_reconst[i] = efun_reconst_i;
    
    arma::vec muO_i(locO.size());
    //argvalsO_i;
    for(arma::uword j = 0; j < locO.size(); ++j){
      muO_i[j] = mu[locO[j]];
    }
    muO[i] = muO_i; 
    evaluesOO[i] = evaluesO_i;
    efunctionsO[i] = efunctionsO_i;
  }
  
  return std::make_tuple(muO, scoresO, CE_scoresO, efunctions, efunctionsO,
                        efun_reconst, evalues, evaluesOO, obs_argvalsO, locOO,sigma2, cov_est);
                      
}


int gcvKneipLiebl(const NumericVector& mu, const std::pair<std::vector<double>, NumericMatrix>& Y_preprocessed, 
                  const std::vector<double>& argvalsO_i, const arma::uvec& locO,
                  const arma::mat& cov_est, double sigma2, const std::string& method, double pev)
{
  NumericMatrix Y(clone(Y_preprocessed.second));
  std::vector<int> complete_rows;
  complete_rows.reserve(Y.nrow());
    
  for (int i = 0; i < Y.nrow(); ++i) {
      if (!NumericVector::is_na(Y(i, Y.ncol() - 1))) {
          complete_rows.push_back(i);
      }
  }

  complete_rows.shrink_to_fit();
    
  NumericMatrix Y_c(complete_rows.size(), Y.ncol());//matrice delle righe complete
  for (size_t i = 0; i < complete_rows.size(); ++i) {
    Y_c.row(i) = Y.row(complete_rows[i]);
  }

  std::vector<double> argvals(Y_preprocessed.first);
  std::vector<size_t> locM;
  locM.reserve(argvals.size());
  for(size_t i = 0; i < argvals.size(); ++i){
    if(arma::find(locO == i, 1).is_empty())
      locM.push_back(i);
  }
  locM.shrink_to_fit();

  NumericMatrix Y_pred(clone(Y_c));
  for(int i = 0; i < Y_pred.nrow(); ++i)
  {
    for(auto& col: locM)
    {
      Y_pred(i,col) = NA_REAL; //mask values to predict
    }
  }
  bool too_few = false;
  //select subset for the presmoothing. observations >= 5 <- secondo me questo in R era inutile
  if(Y_pred.ncol() - locM.size() < 5){Rcout<<"not enough observations for the presmoothing"<<std::endl; too_few = true;}

  std::vector<double> w = quadWeights(argvalsO_i);//giusto
  arma::vec w_arma = arma::conv_to<arma::vec>::from(w);
  arma::mat Wsqrt = arma::diagmat(arma::sqrt(w_arma));//giusto
  arma::mat Winvsqrt = arma::diagmat(1 / arma::sqrt(w_arma));//giusto
  // CovOO
  arma::mat VO = Wsqrt * cov_est.submat(locO,locO) * Wsqrt;//giusto
  arma::vec evaluesO;  //in realtà evaluesO ed eigenvectorsO dovrei averceli già nei data member, devo passarli a gcv
  arma::mat eigenvectorsO;
  arma::eig_sym(evaluesO, eigenvectorsO, VO);//se li passo posso saltare questa parte
  arma::uvec sort_indices = arma::sort_index(evaluesO, "desc");
  evaluesO = evaluesO(sort_indices);
  eigenvectorsO = eigenvectorsO.cols(sort_indices);
  evaluesO.transform([](double val){ return val <= 0 ? 0.0 : val;});//quelli molto piccoli sono sbagliati ... ma molto vicini allo 0
  size_t npcO = arma::sum(evaluesO > 0);
  if(npcO == 0){stop("no eigenvalues greater than 0");}else{npcO--;};//vedi se togliere 1, ricontrolla anche se in eigen ha senso
  double sum_pos_O = arma::sum(evaluesO(arma::find(evaluesO > 0)));
  if (!NumericVector::is_na(pev)) { //riguardare default di pev, in R è NULL, qua sto facendo come se fosse NA_REAL
    size_t npc_index_i = 0;

    arma::vec cumsum_pos_O = arma::cumsum(evaluesO(arma::find(evaluesO > 0)));
    double threshold = pev * sum_pos_O;
    arma::uvec indices = arma::find(cumsum_pos_O >= threshold);

    if(!indices.is_empty())//in realtà non dovrebbe essere mai empty
      {npc_index_i = indices[0]; npcO = npc_index_i;}
    
  }//forse in realtà è giusto rifare questi calcoli da capo perchè sto usando Y_pred??
  //npcO giusto (uno in meno di R)
  arma::mat efunctionsO_i = Winvsqrt * eigenvectorsO.cols(0, npcO);//giusto
  arma::vec evaluesO_i = evaluesO.subvec(0, npcO);//giusto
  arma::mat D_inv = arma::diagmat(1/evaluesO_i);
  arma::mat Z = efunctionsO_i;

  NumericMatrix efun_reconst_i(argvals.size(), npcO + 1);//+1 perchè era un indice
  efun_reconst_i.fill(NA_REAL);    
  arma::mat rows = cov_est.rows(locO);
  for(int k = 0; k < efun_reconst_i.ncol(); ++k)
  {
    for(int r = 0; r < efun_reconst_i.nrow(); ++r){
        arma::vec c = rows.col(r);
        double integral = trapezioidal_rule(efunctionsO_i.col(k) % c, arma::conv_to<arma::vec>::from(argvalsO_i));
        efun_reconst_i(r,k) = integral/evaluesO_i[k];
    }
  }//giusto

  if(too_few)
  {
    Rcout<<"too few. Threshold 0.9"<<std::endl;
    size_t npc_index_i;
    arma::vec cumsum_pos_O = arma::cumsum(evaluesO(arma::find(evaluesO > 0)));
    double threshold = 0.9 * sum_pos_O;
    arma::uvec indices = arma::find(cumsum_pos_O >= threshold);

    if(!indices.is_empty())//in realtà non dovrebbe essere mai empty
      {npc_index_i = indices[0]; npcO = npc_index_i;return npcO;}//npcO = K.pve in R  
  }
  //effettivamente dovrei crearmi la variabile n_comp e usarla al posto di .nrow()

  NumericMatrix rss_mat(Y_pred.nrow(),npcO+1);
  rss_mat.fill(NA_REAL);//se sono arrivata fino a qua c'erano righe complete
  for(int i = 0; i < rss_mat.nrow(); ++i){
    NumericVector Y_cent = Y_pred(i,_) - mu;
    arma::uvec obs_locO(Y_pred.ncol() - sum(is_na(Y_cent)));//giusto
    for(arma::uword j = 0; j < obs_locO.size(); ++j){
      obs_locO[j] = j;
    }
    arma::vec y_cent_arma(obs_locO.n_elem);
    for(const auto& u:obs_locO)
    {
      y_cent_arma[u] = Y_cent[u];
    }
    NumericVector fragmO_presmooth;//rimane empty se non entra in KLAl4
    if(method == "KLAl")//else è KLAl5 e lo salta
    {
      NumericVector y;
      NumericVector x = wrap(argvalsO_i); //argvalsO_i[obs_locO] = argvalsO_i per costruzione!
      for(const auto&u: obs_locO)
      {
        y.push_back(Y_pred(i,u));//stats::na.omit((Y.pred[i,]))
      }
      
      Environment stats = Environment::namespace_env("stats");
      Function smooth_spline = stats["smooth.spline"];
      //smooth splines
      List smooth_fit = smooth_spline(_["y"] = y, _["x"] = x);//S3 object
      Function predict = stats["predict"];
      List result = predict(smooth_fit, _["x"] = x);
      fragmO_presmooth = result["y"];
    }
    //if method == "KLAl4" || method == "KLAl5"
    arma::mat Zcur, ZtZ_sD_inv; arma::vec CE_scoresO_i;
    if(sigma2 == 0.0){sigma2 = 1e-6;}
    if(obs_locO.size() < npcO + 1){
      size_t npcO_new = obs_locO.size(); //ora è una lunghezza, non indice
      Zcur = Z.submat(arma::span(obs_locO[0], obs_locO[obs_locO.size()-1]),arma::span(0, npcO_new - 1));//controlla
      ZtZ_sD_inv = arma::inv(Zcur.t() * Zcur + sigma2 * D_inv.submat(arma::span(0, npcO_new - 1),arma::span(0, npcO_new - 1)));
      CE_scoresO_i = ZtZ_sD_inv * Zcur.t() * y_cent_arma;//arma::vec sono vettori colonna, controlla dimensioni
      size_t new_size = CE_scoresO_i.size() + npcO + 1 - npcO_new;
      CE_scoresO_i.resize(new_size);

      for (size_t i = CE_scoresO_i.size() - (npcO + 1 - npcO_new); i < CE_scoresO_i.size(); ++i) {
          CE_scoresO_i[i] = 0;
      }
    }else{
      Zcur = Z.rows(obs_locO);//giusto a meno di segni
      ZtZ_sD_inv = arma::inv(Zcur.t() * Zcur + sigma2 * D_inv);//giusto a meno di segni
      CE_scoresO_i = ZtZ_sD_inv * Zcur.t() * y_cent_arma;//giusto a meno di segni
    }
    for(size_t k = 1; k <= npcO + 1; ++k)//perchè seq_len include l'estremo destro
    {
      std::tuple<arma::vec, arma::vec> result_tmp = reconstKL_fun(mu, argvals, locO, 
                                      CE_scoresO_i, efun_reconst_i, fragmO_presmooth, k);
      arma::vec y_reconst = std::get<0>(result_tmp);
      double sum = 0.0;
      for(const auto&m : locM)
      {
        sum += (y_reconst[m] - Y_c(i,m))*(y_reconst[m] - Y_c(i,m));//std::pow è inefficiente, cambiarlo anche da altre parti
      }
      rss_mat(i,k-1) = sum;
    }

  }//end for
  std::vector<double> gcv_k_vec; gcv_k_vec.reserve(npcO+1);
  for(size_t i = 0; i < npcO + 1; ++i)
  {
    double sum = 0.0;
    for(int r = 0; r < rss_mat.nrow(); ++r)
    {
      sum += rss_mat(r,i);//colsum
    }
    double denom = (1-(i+1)/static_cast<double>(Y_pred.nrow())) * (1-(i+1)/static_cast<double>(Y_pred.nrow()));//ncompl = Y_pred.nrow()
    gcv_k_vec.push_back(sum/denom);//giusto
  }
  auto it_min = std::min_element(gcv_k_vec.begin(), gcv_k_vec.end());
  auto index_min = std::distance(gcv_k_vec.begin(), it_min);
  // K = index_min + 1, perchè K in R parte da 1, non può essere 0
  return index_min + 1;
}