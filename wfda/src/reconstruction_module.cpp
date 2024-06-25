#include "reco_factory.h"
RCPP_MODULE(reconstruction) {
Rcpp::class_<ReconstructionBase>("ReconstructionBase")
.factory<const std::string&,const NumericMatrix&>(reconstructionFactory)//expose the factory function
.method("reconstruct", &ReconstructionBase::reconstructCurve)
.method("reconst_fcts", &ReconstructionBase::reconst_fcts)
.method("mean",&ReconstructionBase::mean)
.method("cov",&ReconstructionBase::cov);
}



