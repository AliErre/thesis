#include "genetic.h"
RCPP_MODULE(Genetic){
Rcpp::class_<Genetic>("Genetic")
.constructor<List, List, arma::mat, List, arma::vec, arma::uvec>()
.method("multistart", &Genetic::multistart)
.method("get_Pin", &Genetic::get_Pin)
.method("get_vin", &Genetic::get_Vin)
.method("get_P", &Genetic::get_P)
.method("get_v", &Genetic::get_V)
.method("get_blist", &Genetic::get_blist);
}