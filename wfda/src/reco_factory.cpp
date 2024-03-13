# include "reco_factory.h"

//TODO: capire se passare tpoints come const reference o no!!
std::unique_ptr<ReconstructionBase> reconstructionFactory(const std::string& id, const NumericMatrix& Y) {
    if (id == "Kraus") {
        return std::make_unique<ReconstructionKraus>(Y);
    } else if (id == "KLAl") {
        return std::make_unique<ReconstructionKLAl>(Y);
    } else {
        return nullptr
    }
}
// usage example: recon_object <- reconstructionFactory("Kraus", matrix_data, 0.5, 10, t_points_vector)