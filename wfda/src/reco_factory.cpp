# include "reco_factory.h"

//.factory function typically expects a function pointer that returns a raw pointer (T*) rather than a smart pointer like std::unique_ptr. 
//R's memory management model doesn't directly align with C++ smart pointers
ReconstructionBase* reconstructionFactory(const std::string& id, const NumericMatrix& Y) {
    if (id == "Kraus") {
        return new ReconstructionKraus(Y);//must use raw pointers to be compatible to R
    } /*else if (id == "KLAl") {
        return new ReconstructionKLAl(Y);
    }*/ else {
        return nullptr;
    }
}
// usage example: recon_object <- reconstructionFactory("Kraus", matrix_data, 0.5, 10, t_points_vector)