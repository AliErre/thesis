#ifndef GENETIC_HELPER_FUNC_H
#define GENETIC_HELPER_FUNC_H
#include "RcppArmadillo.h"
#include <tuple>
using namespace Rcpp;
bool is_fdPar(const List&);
void my_fRegressArgCheck(List&, List&, List&, const Function&, const Function&, const Function& fdPar);
List predict_fRegress(const std::tuple<List, List, List, List, List, arma::mat, arma::vec>&, const List&, const arma::vec&,
                      const Function&, const Function&);
std::tuple<List, List, List, List, List, arma::mat, arma::vec> 
weighted_fRegress(List&, List&, List&, const Nullable<List>&, bool, const Function& ,const Function&, 
                  const Function&, const Function&, const Function&, const Function&, const Function&, 
                  const Function&, const Function&, const Function&);

#endif //GENETIC_HELPER_FUNC_H