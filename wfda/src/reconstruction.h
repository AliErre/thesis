#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H
#include <memory>
#include <vector>
#include <Rcpp.h>
using namespace Rcpp;

// TODO: return type of reconstructCurves and argument types
//       controllare alpha (Kraus) di che tipo dev'essere e decidere un default
//             //      K   (KLAl)                     //
// in R, reconstKraus_fun ritorna una lista. Cosa ritornare in C++? tupla?
// aggiungere i parametri di KLNoAl

// con Rcpp non esiste variadic templates, se voglio usare la factory non posso passare diversi numeri di parametri
// quando derived class non ha quel parametro (ex: Kraus non ha K) passarlo come NA
// implementare gestione di NA, forse i tipi di Rcpp li accettano senza problemi
class ReconstructionBase {
    public:
        ReconstructionBase(const NumericMatrix& Y, const Numeric alpha = 0.0, 
                           const Integer K = 0, const NumericVector& t_points = NumericVector(0),
                           Integer nRegGrid = 0, Integer maxBins = 0) 
                          : m_Y(Y), m_alpha(alpha), m_K(K), m_t_points(t_points), m_nRegGrid(nRegGrid), m_maxBins(maxBins)
                          { m_reconst_fcts = find_obs_inc(Y);}  //constructor

        // specialize
        virtual ~ReconstructionBase() = default; // polymorphism => need virtual destructor
        virtual void reconstructCurves() = 0; // override in derived classes
        virtual Numeric alpha() const = 0;
        virtual Integer K() const = 0;
        virtual NumericVector t_points() const = 0;
        virtual Integer max_bins() const = 0;
        virtual Integer n_reg_grid() const = 0;

        // same for all derived
        std::vector<bool> find_obs_inc(const NumericMatrix& Y) const; //Rcpp will deal with return type
        std::vector<bool> reconst_fcts() {return m_reconst_fcts;} //getter

    protected:
        NumericMatrix m_Y; // curves.train in R
        std::vector<size_t> m_reconst_fcts; //vector of indices of curves to reconstruct 
        Numeric m_alpha;
        Integer m_K; //check type
        NumericVector m_t_points; //check type
        Integer m_maxBins;
        Integer m_nRegGrid;
};


class ReconstructionKraus : public ReconstructionBase{
    public:
        //constructor for Kraus
        ReconstructionKraus(const NumericMatrix& Y, const Numeric alpha, const Integer K, 
                            const NumericVector& t_points, Integer nRegGrid = 0, Integer maxBins = 0) : 
                            ReconstructionBase(Y, alpha, K, t_points, nRegGrid, maxBins) //per unique_ptr
        void reconstructCurves() override;//metodo che sarà chiamato da R

        std::vector<double> meanKraus() const;
        std::vector<std::vector<double>> covKraus() const;
        

        // getter e setter mean e cov
        std::vector<double> mean() const { return m_mean; } //returns value computed by meanKraus
        std::vector<std::vector<double>> cov() const { return m_cov; } //returns value computed by covKraus

        double gcvKraus(const std::vector<std::vector<double>>& covMat, const std::vector<double>& meanVec, 
                        const NumericMatrix& X, const bool M_bool_vec, const Numeric& alpha) const;
        TUPLA reconstKraus_fun(const std::vector<std::vector<double>>& covMat, 
                               const std::vector<double>& X_cent_vec, 
                               const Numeric& alpha); //ritorna una lista => how to deal?

        //override
        Numeric alpha() const override {return m_alpha; }
        Integer K() const override { Rcout << "no K for Kraus method" << std::endl; }
        NumericVector t_points() const override { Rcout << "no t_points for Kraus method" << std::endl; }
        Integer max_bins() const override { Rcout << "no bins for Kraus" << std::endl; }
        Integer n_reg_grid() const override { Rcout << "no n_reg_grid for Kraus" << std::endl; }

    private:
        std::vector<double> m_mean; // set to NA default? capire se mettere return type NumericVector
        std::vector<std::vector<double>> m_cov; //default? return type NumericMatrix?

};

class ReconstructionKLAl : public ReconstructionBase{
    public: //const Integer K& = 0, const NumericVector& t_points = NumericVector(0) <- nel constructor
        ReconstructionKLAl(const NumericMatrix& Y, const Numeric alpha, const Integer K, 
                           const NumericVector& t_points, Integer nRegGrid = 0, Integer maxBins = 0) : 
                           ReconstructionBase(Y, alpha, K, t_points, nRegGrid, maxBins) //per unique_ptr
        void reconstructCurves() override;

        Integer K() const override { return m_K; }
        NumericVector t_points() const override { return m_t_points; } //controllare se override va messo qua?
        Integer max_bins() const override { return m_maxBins; }
        Integer n_reg_grid() const override { return m_nRegGrid; }
        Integer n_reg_grid() const override { Rcout << "no n_reg_grid for Kraus" << std::endl; }
        Numeric alpha() const override { Rcout<< "no alpha for KL method" << std::endl; }

    private:

};

#endif RECONSTRUCTION_H