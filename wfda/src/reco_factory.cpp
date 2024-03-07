# include "reco_factory.h"

//TODO: ai costruttori vanno passati gli argomenti
std::unique_ptr<ReconstructionBase> reconstructionFactory(const std::string& id, const NumericMatrix& Y, 
                                                          Numeric alpha, Integer K = 0, const NumericVector& tpoints = NumericVector(0),
                                                          Integer nRegGrid = 0, Integer maxBins = 0) {
    if (id == "Kraus") {
        return std::make_unique<ReconstructionKraus>(Y, alpha, K, tpoints);
    } else if (id == "KLAl") {
        return std::make_unique<ReconstructionKLAl>(Y, alpha, K, tpoints);
    } else {
        return nullptr
    }
}
// usage example: recon_object <- reconstructionFactory("Kraus", matrix_data, 0.5, 10, t_points_vector)