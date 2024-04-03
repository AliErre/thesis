# include "reco_factory.h"

//.factory function typically expects a function pointer that returns a raw pointer (T*) rather than a smart pointer like std::unique_ptr. 
//R's memory management model doesn't directly align with C++ smart pointers
ReconstructionBase* reconstructionFactory(const std::string& id, const NumericMatrix& Y) {
    if (id == "Kraus") {
        return new ReconstructionKraus(id, Y);//must use raw pointers to be compatible to R
    } else if (id == "KLAl" || id == "KLNoAl") {
        return new ReconstructionKL(id, Y);
    } else if (id == "Extrapolation"){
        return new ReconstructionExtrapolation(id, Y);//passare al costruttore anche id e risolvo tutti i prob
    }
    else {
        return nullptr;
    }
}
// usage example: recon_object <- reconstructionFactory("Kraus", matrix_data, 0.5, 10, t_points_vector)