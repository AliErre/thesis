#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat getEigenValues(arma::mat M) {
    arma::vec evalues;
    arma::mat evectors;
    arma::eig_sym(evalues, evectors,M);
    return evectors;
}
// [[Rcpp::export]]
List try_list(const List& r)
{
  List new_list = clone(r);
  for(int i = 0 ; i < r.length(); i++)
  {
    int elem = new_list[i];
    new_list[i] = elem + 1;
  }
  return r;
}

// [[Rcpp::export]]
NumericMatrix test_fda(List& xxfdjk, List& betabasisj, List& betabasisk, NumericVector& range){
  Environment fda = Environment::namespace_env("fda");
  Function sum_fd = fda["sum.fd"];
  Function inprod = fda["inprod"];
  List wtfdjk = sum_fd(xxfdjk);
  NumericMatrix Cmatjk = inprod(_["fdobj1"] = betabasisj, _["fdobj2"] = betabasisk, _["rng"] = range, _["wtfd"] = wtfdjk);
  return Cmatjk;
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
RObject sum_cpp(arma::vec vec1) {
  Function sum = Environment::base_env()["sum"];
  return sum(vec1);
}
/*// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::uvec test_set(arma::vec vec1, arma::vec vec2)
{
  arma::uvec test;
  for(const auto& t:vec2)
  {
      arma::uvec indices = arma::find(vec1 == t);
      test.insert_rows(test.n_elem, indices);
  }
  return test;
}
// [[Rcpp::export]]
arma::vec sample_cpp(const arma::vec& vec)
{
    //RNGkind(sample.kind = "Rounding")
    arma::vec ret = RcppArmadillo::sample(vec, 3, false);
    return ret;
}*/
/*// [[Rcpp::export]]*/
arma::vec eigenvalues(arma::mat V) {
    // Compute eigenvalues of a symmetric matrix
    arma::vec eigvals = arma::eig_sym(V);

    return eigvals;
}

/*// [[Rcpp::export]]*/
NumericMatrix inprod_cpp(List& fd1, List& fd2, NumericVector& range, List& w)
{
  Environment fda = Environment::namespace_env("fda");
  Function inprod = fda["inprod"];
  Function int2Lfd = fda["int2Lfd"];
  NumericMatrix ret = inprod(_["fdobj1"] = fd1, _["fdobj2"] = fd2, _["Lfdobj1"] = int2Lfd(0), _["Lfdobj2"] = int2Lfd(0), _["rng"] = range, _["wtfd"] = w);
  return ret;
}

/*// [[Rcpp::export]]*/
List vec(const NumericMatrix& G0, const NumericVector& row_vec, const NumericVector& col_vec, const NumericVector& weights)
{
  Environment mgcv = Environment::namespace_env("mgcv");
  Function gam = mgcv["gam"];
  Function predict_gam = mgcv["predict.gam"];
  NumericVector g0(G0.begin(),G0.end());
  int nbasis = 10;
  std::string formula = "G0 ~ te(row_vec, col_vec, k = " + std::to_string(nbasis) + ")";
  Formula f = Formula(formula);
  DataFrame data = DataFrame::create(_["G0"] = g0, _["row_vec"] = row_vec, _["col_vec"] = col_vec);
  List gamModel = gam(_["formula"] = f, _["data"] = data, _["weights"] = weights);//fit gam model
  Rcout<<"fitted gam model in smooth_cov"<<std::endl;
  DataFrame newdata = DataFrame::create(Named("row_vec") = row_vec, Named("col_vec") = col_vec);//data for prediction

  NumericVector predictions = predict_gam(gamModel, _["newdata"] = newdata);

  return gamModel;
}

/*// [[Rcpp::export]]*/
List compute(Formula& f, NumericVector& Y, NumericVector& d_vec)
{   
    Environment stats = Environment::namespace_env("stats");
    Function gam = stats["gam"];
    DataFrame data = DataFrame::create(Named("Y") = Y, Named("d.vec") = d_vec);
    List result = gam(f, Named("data") = data);
    return result;
}
// [[Rcpp::export]]
NumericVector check(NumericVector& Y_first, int d){
  NumericVector row_vec;
  for(const auto& y:Y_first)
  {
    for(int j = 0; j < d; ++j){
      row_vec.push_back(y);}
  }
  return row_vec;
}
/*
// [[Rcpp::export]]
std::vector<size_t> generateIdVec(const std::vector<std::vector<double>>& Ly) {
    for(const auto& L:Ly)
    {
      for(const auto& l:L)
        Rcout<<l<<"\t";
      Rcout<<std::endl;
    }
    std::vector<size_t> id_vec;
    size_t total_size = 0;

    // Calculate the total size needed for id_vec
    for (size_t i = 0; i < Ly.size(); ++i) {
        total_size += Ly[i].size();
    }

    // Reserve memory for id_vec
    id_vec.reserve(total_size);

    for (size_t i = 0; i < Ly.size(); ++i) {
        Rcout<<"i:"<<i<<"\t";
        for (size_t j = 0; j < Ly[i].size(); ++j) {
            id_vec.push_back(i); // Assigning identifiers based on the index of the outer vector
            Rcout<<i<<"\t";
        }
    }
    for(const auto&id: id_vec)
    {
      Rcout<<id<<"\t";
    }
    return id_vec;
}*/