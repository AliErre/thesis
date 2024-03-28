#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
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