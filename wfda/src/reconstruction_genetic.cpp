#include "genetic.h"
RCPP_MODULE(Genetic){
Rcpp::class_<Genetic>("Genetic")
.constructor<List, List, arma::mat, List, arma::vec, arma::uvec>()
.method("multistart", &Genetic::multistart);
}