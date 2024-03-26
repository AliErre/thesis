#include <RcppArmadillo.h>
#include "helper_functions.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::export]]
List reconstKL_fun(const NumericVector& mu, const std::vector<double>& argvals, const arma::uvec& locO, 
                   const arma::vec& scoresO, const NumericMatrix& efunc_r, const NumericVector& fragmO, int k)
{
  //non so se esportarla, o se la esporto devo avere a disposzione tutti sti dati comunque
  int K = std::min(k, efunc_r.ncol());
  arma::mat efunc_r_arma = as<arma::mat>(efunc_r);
  arma::vec mu_arma = as<arma::vec>(mu);
  arma::vec reconstr = mu_arma.t() + (arma::mat(scoresO)).rows(0,K-1).t() * (efunc_r_arma.cols(0,K-1)).t();//metterlo in arma::mat dovrebbe permettermi di avere un vettore riga
  if(fragmO.size() != 0)//metodo4
  {
    //arma::uword min = locO.min();//locO min sarà sempre 0?
    //reconstr.subvec(0, min) += fragmO[0] - reconstr[min];
    arma::uword max = locO.max();
    reconstr.subvec(max + 1, argvals.size() - 1) += fragmO[fragmO.size() - 1] - reconstr[max];
    reconstr(locO) = as<arma::vec>(fragmO);//controlla poi che anche a me fragmO venga di dimensione locO.size()
  }

  //weights null, metterlo qua? in R ritorna un NULL
  return List::create(Named("y_reconst") = reconstr);//ritorna anche argvals ma mi sembra una cosa scema visto che è l'argomento con cui è chiamata
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
                                             bool binning = false, int max_bins){ //vedi se aggiungere il default = 1000 per maxBins

  //crea una copia
  std::vector<std::tuple<int, double, double>> y_data_complete(find_complete_tuple(y_data));
  std::set<double> bins; 

  std::vector<int> ids;//controllare se qua gli elementi partono da 0 o da 1. Per ora do per scontato partano da 0
  ids.reserve(y_data_complete.size());
  for(const auto& t: y_data_complete)
  {
    ids.push_back(std::get<0>(t));//newid. deve essere il vettore di indici di riga di Y
    bins.insert(std::get<1>(t));
  }//bins è per forza uguale a t_points
  //ids sarà già ordinato, non c'è bisogno di chiamare sort
  int nobs = std::distance(ids.begin(),std::unique(ids.begin(),ids.end()));//corretto
  bool condition = binning && bins.size() > max_bins;//false
  std::vector<double> binvalues = make_bins(bins, max_bins,condition);//bins sarà modificato, chiamo make_bins con la reference
  //associa a std::get<1>(t) un bin => crea classi di .index

  //estrai vettore dal secondo elemento della tupla
  std::vector<double> index_vector;
  index_vector.reserve(y_data_complete.size());
  for(const auto& t:y_data_complete)
  {
    index_vector.push_back(std::get<1>(t));
  }
  index_vector.shrink_to_fit(); //non dovrebbe servire perchè la dimensione è quella precisa

  std::vector<std::pair<double, double>> classes = cut(index_vector, bins);//this is newindex in R
  //associate each class to a column index
  //index classes by the ranking of their lower bound
  /*std::vector<std::pair<double,double>> sorted_classes(classes.begin(), classes.end());
  std::sort(sorted_classes.begin(),sorted_classes.end(),[](const std::pair<double, double>& a, const std::pair<double, double>& b)
                                                        { return a.first < b.first;});
                                                        //lambda specifies criterium to order: order on lower bound
  */
  std::set<std::pair<double,double>> sorted_classes(classes.begin(),classes.end());
  //associate index
  std::vector<int> column_indices(classes.size()); //vedi se trasformare in size_t       

  for (int i = 0; i < classes.size(); ++i) {
      auto it = std::find(sorted_classes.begin(), sorted_classes.end(), classes[i]);//lo trova per forza. find su un vettore restituisce iteratore
      column_indices[i] = std::distance(sorted_classes.begin(), it);
  }

  //ora costruisci Y(ids[i],column_index[i]) = y_data_complete[i];
  //controlla di aver costruito ids e column_index a partire da y_data_complete
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


std::pair<NumericMatrix,NumericVector> smooth_cov(const NumericMatrix& Y_second,const NumericMatrix& Y_tilde,const NumericVector& Y_first,
                         int d, int i, int nbasis)
{
  NumericMatrix cov_mean(d, d);
  NumericMatrix cov_count(d,d);
  std::fill(cov_count.begin(),cov_count.end(),0.0);
  NumericMatrix cov_sum(d,d);
  std::fill(cov_sum.begin(),cov_sum.end(),0.0);


  // Calculate cov.mean
  for (int idx = 0; idx < i; ++idx) {

    IntegerVector obs_points = seq_len(d) - 1;//0 based index
    IntegerVector obs_points = obs_points[!is_na(Y_second(idx, _))];

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

  // Calculate G.0
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


  // Calculate npc.0
  NumericVector row_vec(d * Y_first.size());
  for(int i = 0; i < Y_first.size();++i)
  {
    for(int j = 0; j < d; ++j)
    row_vec[j] = Y_first[i];
  }


  NumericVector col_vec(Y_first.size() * d);
  for(int i = 0; i < d; ++i) {
    std::copy(Y_first.begin(), Y_first.end(), col_vec.begin() + i * Y_first.size());
  }

  Environment mgcv = Environment::namespace_env("mgcv");
  Function gam = mgcv["gam"];
  Function predict_gam = mgcv["predict.gam"];
  NumericVector g0(G0.begin(),G0.end());
  NumericVector weights(cov_count.begin(),cov_count.end());

  DataFrame gamData = DataFrame::create(Named("G0") = g0, 
                                        Named("row_vec") = row_vec, 
                                        Named("col_vec") = col_vec, 
                                        Named("weights") = weights);
  
  List gamModel = gam(Formula("G0 ~ te(row_vec, col_vec, k = " + std::to_string(nbasis) + ")"), gamData, _["weights"] = weights);//fit gam model
  DataFrame newdata = DataFrame::create(Named("row.vec") = row_vec, Named("col.vec") = col_vec);//data for prediction
  NumericVector predictions = predict_gam(gamModel, _["newdata"] = newdata);

  //reshape predictions vector into matrix
  NumericMatrix npc0(d, d);
  std::copy(predictions.begin(), predictions.end(), npc0.begin());
  NumericMatrix npc0_sym(d, d);

  // symmetrization
  for (int i = 0; i < d; ++i) {
    for (int j = 0; j < d; ++j) {
          npc0_sym(i, j) = (npc0(i, j) + npc0(j, i)) / 2.0;
    }
  }

  return std::make_pair(npc0_sym,diag_G0);
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
        ret[0] = 0.0; // Assuming the first weight for midpoint is 0 as in the R example
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
  for(size_t i = 0; i < v.size(); ++i)
  {
    sum += v[i]*w[i];
    sum_w += w[i];//this does not necessarily sum to 1 ?
  }
  return sum/sum_w; 
}


double trapezioidal_rule(const arma::vec& values, const arma::vec& knots){
  double integral = 0.0;
  int n = knots.size();

  for (arma::uword i = 1; i < n; ++i) {
      integral += (knots[i] - knots[i - 1]) * (values[i] + values[i - 1]) / 2.0;
  }

  return integral;
}

std::tuple<List, List, List, arma::mat, List, List, arma::vec, List, List, List, double, arma::mat> eigen(const std::vector<double>& argvals, const std::vector<std::vector<double>>& argvalsO,const NumericMatrix& npc_0, 
           double pev, const NumericMatrix& Y_pred, 
           const NumericVector& mu, const NumericVector& diagG0, bool CEScores, const IntegerVector& reconst_fcts) {
  
  // numerical integration for calculation of eigenvalues
  std::vector<double> w = quadWeights(argvals);//default is trapezioidal
  arma::vec w_arma = arma::conv_to<arma::vec>::from(w);
  arma::mat Wsqrt = arma::diagmat(arma::sqrt(w_arma));
  arma::mat Winvsqrt = arma::diagmat(1 / arma::sqrt(w_arma));
  arma::mat V = Wsqrt * as<arma::mat>(npc_0) * Wsqrt;

  arma::vec evalues;
  arma::mat eigenvectors;
  arma::eig_sym(evalues, eigenvectors, V);
  evalues.transform([](double val) { return val <= 0 ? 0.0 : val; });//ora tutti gli autovalori sono >= 0
  int npc = arma::sum(evalues > 0) - 1;//conta valori positivi -> poi lo devo usare come indice di posizione, quindi tolgo 1
  double sum_pos = arma::sum(evalues(arma::find(evalues > 0)));

  if (!NumericVector::is_na(pev)) { //riguardare default di pev, in R è NULL, qua sto facendo come se fosse NA_REAL
    arma::uword npc_index;

    if (sum_pos > 0) {
      arma::vec cumsum_pos = arma::cumsum(evalues(arma::find(evalues > 0))); // Cumulative sum of positive eigenvalues
      double threshold = pev * sum_pos; // Calculate the threshold value
      arma::uvec indices = arma::find(cumsum_pos >= threshold);

      if(!indices.is_empty())//in realtà non dovrebbe essere mai empty
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

  double lower_bound = argvals[0] + 0.25 * T_len;//can use it only on an ordered container, so argvals must be order (which it is)
  auto it = std::lower_bound(argvals.begin(), argvals.end(), lower_bound);
  int T1_min = std::distance(argvals.begin(), it); //in R which returns an index

  double upper_bound = argvals[argvals.size()] - 0.25 * T_len;
  auto it = std::upper_bound(argvals.begin(), argvals.end(), upper_bound);
  if(it != argvals.begin())
  {
    --it;
  }

  int T1_max = std::distance(argvals.begin(), it);// = 0 se sono tutti maggiori di upper_bound. Ma non dovrebbe succedee
  arma::vec diag_diff = as<arma::vec>(diagG0) - arma::diagvec(cov_est);
  arma::vec sub_diag = diag_diff.rows(T1_min, T1_max);
  std::vector<double> argvals_subset(argvals.begin() + T1_min, argvals.begin() + T1_max);
  std::vector<double> w2 = quadWeights(argvals_subset);//default: trapezioidal

  //sigma
  double sigma2 = std::max(weighted_mean(argvals_subset, w2), 0.0);//lei aveva messo na.rm ma secondo me non dovrebbero esserci NA
  
  //computations for observed fragments
  //cambia queste definizioni NOTA BENE NOTA BENE NOTA BENE?
  List muO(reconst_fcts.size()), scoresO(reconst_fcts.size()), CE_scoresO(reconst_fcts.size()),
           evaluesOO(reconst_fcts.size()), efunctionsO(reconst_fcts.size()), efun_reconst(reconst_fcts.size()),
           obs_argvalsO(reconst_fcts.size()), locOO(reconst_fcts.size());
  
  
  for (int i = 0; i < reconst_fcts.size(); ++i) {//argvalsO should be of size reconst_fcts (observed_period in myfpca)
    // Numerical integration for calculation of eigenvalues
    std::vector<double> w = quadWeights(argvalsO[i]);
    arma::vec w_arma = arma::conv_to<arma::vec>::from(w);
    arma::mat Wsqrt = arma::diagmat(arma::sqrt(w_arma));
    arma::mat Winvsqrt = arma::diagmat(1 / arma::sqrt(w_arma));

    arma::uvec locO(argvalsO[i].size());//0-indexed. su R è 1-indexed
    for(arma::uword j = 0; j < locO.size(); ++j)
    {
      locO[j] = j; //è come match() su R, per la natura di argvalsO e argvals
    }
    locOO[i] = locO;

    // CovOO
    arma::mat VO = Wsqrt * cov_est.submat(locO,locO) * Wsqrt;
    arma::vec evaluesO;
    arma::mat eigenvectorsO;
    arma::eig_sym(evaluesO, eigenvectorsO, VO);
    evaluesO.transform([](double val){ return val <= 0 ? 0.0 : val;});
    int npcO = arma::sum(evaluesO > 0) - 1;

    double sum_pos_O = arma::sum(evaluesO(arma::find(evaluesO > 0)));

    if (!NumericVector::is_na(pev)) { //riguardare default di pev, in R è NULL, qua sto facendo come se fosse NA_REAL
      arma::uword npc_index_i;

      if (sum_pos > 0) {
        arma::vec cumsum_pos_O = arma::cumsum(evaluesO(arma::find(evaluesO > 0))); // Cumulative sum of positive eigenvalues
        double threshold = pev * sum_pos_O; // Calculate the threshold value
        arma::uvec indices = arma::find(cumsum_pos_O >= threshold);

        if(!indices.is_empty())//in realtà non dovrebbe essere mai empty
          {npc_index_i = indices[0]; npcO = static_cast<int>(npc_index_i);}
      }
    }

    arma::mat efunctionsO_i = Winvsqrt * eigenvectorsO.cols(0, npcO);
    arma::vec evaluesO_i = evaluesO.subvec(0, npcO);

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
      if (obs_locO.size() < npcO + 1) {//+1 perchè era un indice e ora lo uso come dimensione
        npcO = obs_locO.size();
        size = true;//se è qua è una size e non un index
      }
      arma::mat Zcur = Z.submat(arma::span(obs_locO[0],obs_locO[obs_locO.size()-1]),arma::span(0, npcO - 1*CEScores*size));//controlla
      arma::mat ZtZ_sD_inv = arma::inv(Zcur.t() * Zcur + sigma2 * D_inv.submat(arma::span(0, npcO - 1*CEScores*size),arma::span(0, npcO - 1*CEScores*size)));
      arma::vec CE_scoresO_i = ZtZ_sD_inv * Zcur.t() * Y_cent_arma(no_na);//arma::vec sono vettori colonna
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
    NumericMatrix efun_reconst_i(argvals.size(), npcO + 1*(!true), NA_REAL);
    for(int k = 0; k < efun_reconst_i.ncol(); ++k)
    {
      NumericVector c = efun_reconst_i(_,k);
      arma::vec c_arma = as<arma::vec>(c);
      double integral = trapezioidal_rule(efunctionsO_i.col(k) % c_arma, arma::conv_to<arma::vec>::from(argvalsO[i]));
      for(int r = 0; r < efun_reconst_i.nrow(); ++r)
       efun_reconst_i(r,k) = integral/evaluesO_i[k];
    }

    efun_reconst[i] = efun_reconst_i;
    
    // Store other results
    arma::vec muO_i(locO.size());
    //argvalsO_i;
    for(arma::uword j = 0; j < locO.size(); ++j){
      muO_i[j] = mu[locO[j]];
    }
    muO[i] = muO_i; 
    evaluesOO[i] = evaluesO_i;
    efunctionsO[i] = efunctionsO_i;
  }
  //mancano cose da ritornare a myfpca
  //vedi se settare data member poi da myfpca a seconda dell uso che ne deve fare KLAl
  return std::make_tuple(muO, scoresO, CE_scoresO, efunctions, efunctionsO,
                        efun_reconst, evalues, evaluesOO, obs_argvalsO, locOO,sigma2, cov_est);/*List::create(Named("muO") = muO,
                      Named("scoresO") = scoresO,
                      Named("CE_scoresO") = CE_scoresO,
                      Named("efunctions") = efunctions,
                      Named("efunctionsO") = efunctionsO,
                      Named("efun_reconst") = efun_reconst,
                      Named("evalues") = evalues,
                      Named("evaluesO") = evaluesOO,
                      Named("sigma2") = sigma2);*/
                      
}



//per ora ritorna un int ma forse K è double, controlla
int gcvKneipLiebl(const NumericVector& mu, const std::pair<std::vector<double>, NumericMatrix>& Y_preprocessed, 
                  const std::vector<double>& argvalsO_i, const NumericVector& muO, const arma::uvec& locO,
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
  std::vector<arma::uword> locM;
  for(arma::uword i = 0; i < argvals.size(); ++i){
    if(arma::find(locO == i, 1).is_empty())
      locM.push_back(i);
  }

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
  //notazione: muO è m_muO[i]

  std::vector<double> w = quadWeights(argvalsO_i);
  arma::vec w_arma = arma::conv_to<arma::vec>::from(w);
  arma::mat Wsqrt = arma::diagmat(arma::sqrt(w_arma));
  arma::mat Winvsqrt = arma::diagmat(1 / arma::sqrt(w_arma));
  // CovOO
  arma::mat VO = Wsqrt * cov_est.submat(locO,locO) * Wsqrt;
  arma::vec evaluesO;  //in realtà evaluesO ed eigenvectorsO dovrei averceli già nei data member, devo passarli a gcv
  arma::mat eigenvectorsO;
  arma::eig_sym(evaluesO, eigenvectorsO, VO);//se li passo posso saltare questa parte
  evaluesO.transform([](double val){ return val <= 0 ? 0.0 : val;});
  int npcO = arma::sum(evaluesO > 0) - 1;//vedi se togliere 1, ricontrolla anche se in eigen ha senso
  double sum_pos_O = arma::sum(evaluesO(arma::find(evaluesO > 0)));
  if (!NumericVector::is_na(pev)) { //riguardare default di pev, in R è NULL, qua sto facendo come se fosse NA_REAL
    arma::uword npc_index_i;

    arma::vec cumsum_pos_O = arma::cumsum(evaluesO(arma::find(evaluesO > 0))); // Cumulative sum of positive eigenvalues
    double threshold = pev * sum_pos_O; // Calculate the threshold value
    arma::uvec indices = arma::find(cumsum_pos_O >= threshold);

    if(!indices.is_empty())//in realtà non dovrebbe essere mai empty
      {npc_index_i = indices[0]; npcO = static_cast<int>(npc_index_i);}
    
  }//forse in realtà è giusto rifare questi calcoli da capo perchè sto usando Y_pred??

  arma::mat efunctionsO_i = Winvsqrt * eigenvectorsO.cols(0, npcO);//dovrebbe esseere giusto aver sottratto 1 prima
  arma::vec evaluesO_i = evaluesO.subvec(0, npcO);
  arma::mat D_inv = arma::diagmat(1/evaluesO_i);
  arma::mat Z = efunctionsO_i;

  //reconstructive eigenfunctions
  NumericMatrix efun_reconst_i(argvals.size(), npcO + 1, NA_REAL);//+1 perchè era un indice
  for(int k = 0; k < efun_reconst_i.ncol(); ++k)
  {
    NumericVector c = efun_reconst_i(_,k);
    arma::vec c_arma = as<arma::vec>(c);
    double integral = trapezioidal_rule(efunctionsO_i.col(k) % c_arma, arma::conv_to<arma::vec>::from(argvalsO_i));
    for(int r = 0; r < efun_reconst_i.nrow(); ++r)
      efun_reconst_i(r,k) = integral/evaluesO_i[k];
  }

  if(too_few)
  {
    arma::uword npc_index_i;
    arma::vec cumsum_pos_O = arma::cumsum(evaluesO(arma::find(evaluesO > 0)));
    double threshold = 0.9 * sum_pos_O; // Calculate the threshold value
    arma::uvec indices = arma::find(cumsum_pos_O >= threshold);

    if(!indices.is_empty())//in realtà non dovrebbe essere mai empty
      {npc_index_i = indices[0]; npcO = static_cast<int>(npc_index_i);return npcO;}//npcO = K.pve in R  
  }
  //effettivamente dovrei crearmi la variabile n_comp e usarla al posto di .nrow()

  NumericMatrix rss_mat(Y_pred.nrow(),npcO+1,NA_REAL);//se sono arrivata fino a qua c'erano righe complete
  for(int i = 0; i < rss_mat.nrow(); ++i){
    NumericVector Y_cent = Y_pred(i,_) - mu;//aggiungi mu alle cose ricevute da gcvKneipLiebl. forse devo creare mu come riga, non so se cosi fa la differenza
    arma::uvec obs_locO(Y_pred.ncol() - sum(is_na(Y_cent)));//giusto
    for(arma::uword j = 0; j < obs_locO.size(); ++j){
      obs_locO[j] = j;
    }
    arma::vec y_cent_arma(obs_locO.size());//controlla che questa inizializzazione abba senso
    for(const auto&u: obs_locO)
    {
      y_cent_arma[u] = (Y_cent[u]-mu[u]);//vedere se a mu puoi accedere con uword
    }
    NumericVector fragmO_presmooth;//rimane empty se non entra in KLAl4
    if(method == "KLAl4")
    {
      NumericVector y;
      
      NumericVector x(argvalsO_i.begin(),argvalsO_i.end()); //argvalsO_i[obs_locO] = argvalsO_i per costruzione!
      for(const auto&u: obs_locO)
      {
        y.push_back(Y_pred(i,u));//stats::na.omit((Y.pred[i,]))
      }
      //DataFrame data = DataFrame::create(Named("x")=x, Named("y") = y);
      Environment stats = Environment::namespace_env("stats");
      Function smooth_spline = stats["smooth.spline"];
      //smooth splines
      List smooth_fit = smooth_spline(Named("y") = y, Named("x") = x);//S3 object
      Function predict = stats["predict"];
      List result = predict(smooth_fit, x);//dubbi
      fragmO_presmooth = result[2];
    }
    arma::mat Zcur, ZtZ_sD_inv; arma::vec CE_scoresO_i, CE_scoresO;
    if(method == "KLAl4" || method == "KLAl5")
    {
      if(sigma2 == 0.0){sigma2 = 1e-6;}
      if(obs_locO.size() < npcO + 1){
        int npcO_new = obs_locO.size(); //ora è una lunghezza, non indice
        Zcur = Z.submat(arma::span(obs_locO[0], obs_locO[obs_locO.size()-1]),arma::span(0, npcO_new - 1));//controlla
        ZtZ_sD_inv = arma::inv(Zcur.t() * Zcur + sigma2 * D_inv.submat(arma::span(0, npcO_new - 1),arma::span(0, npcO_new - 1)));
        CE_scoresO_i = ZtZ_sD_inv * Zcur.t() * y_cent_arma;//arma::vec sono vettori colonna, controlla dimensioni
        CE_scoresO(CE_scoresO_i.size() + npcO + 1 -npcO_new);//perchè npcO è un indice
      }else{
        Zcur = Z.rows(obs_locO);
        ZtZ_sD_inv = arma::inv(Zcur.t() * Zcur + sigma2 * D_inv);
        CE_scoresO_i = ZtZ_sD_inv * Zcur.t() * y_cent_arma;
      }
      for(int k = 1; k <= npcO + 1; ++k)//perchè seq_len include l'estremo destro e parte da 1, io avevo usato npcO come indice quindi tolto 1
      {
        List result_tmp = reconstKL_fun(mu, argvals, locO, 
                                        CE_scoresO, efun_reconst_i, fragmO_presmooth, k);
        arma::vec y_reconst = result_tmp["y_reconst"];
        double sum = 0.0;
        for(const auto&m : locM)
        {
          sum += (y_reconst[m] - Y_c(i,m))*(y_reconst[m] - Y_c(i,m));//std::pow è inefficiente, cambiarlo anche da altre parti
        }
        rss_mat(i,k-1) = sum;
      }

    }
  }//end for
  std::vector<double> gcv_k_vec; gcv_k_vec.reserve(npcO+1);
  for(int i = 0; i < npcO + 1; ++i)
  {
    double sum = 0.0;
    for(int r = 0; r < rss_mat.nrow(); ++r)
    {
      sum += rss_mat(r,i);//colsum
    }
    double denom = (1-(i+1)/Y_pred.nrow()) * (1-(i+1)/Y_pred.nrow());//ncompl = Y_pred.nrow()
    gcv_k_vec.push_back(sum/denom);
  }

  auto it_min = std::min_element(gcv_k_vec.begin(), gcv_k_vec.end());
  auto index_min = std::distance(gcv_k_vec.begin(), it_min);
  // K = index_min + 1, perchè K in R parte da 1, non può essere 0
  return index_min + 1;
}