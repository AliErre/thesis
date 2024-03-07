#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H
#include <memory>
#include <vector>
#include <Rcpp.h>
using namespace Rcpp;

// TODO: return type of reconstructCurves and argument types
//       controllare K (KLAl) di che tipo dev'essere e decidere un default
// in R, reconstKraus_fun ritorna una lista. Cosa ritornare in C++? tupla?
// aggiungere i parametri di KLNoAl

// con Rcpp non esiste variadic templates, se voglio usare la factory non posso passare diversi numeri di parametri
// quando derived class non ha quel parametro (ex: Kraus non ha K) passarlo come NA
// implementare gestione di NA, forse i tipi di Rcpp li accettano senza problemi
class ReconstructionBase {
    public:
    //nota: references point the memory area of the R object!!
        ReconstructionBase(const NumericMatrix& Y, double alpha = 0.0, 
                           int K = 0, NumericVector& t_points = NumericVector(),
                           int nRegGrid = 0, int maxBins = 0) 
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
        std::vector<int> find_obs_inc(const NumericMatrix& Y) const; //Rcpp will deal with return type
        IntegerVector reconst_fcts() const;//getter

    protected:
        NumericMatrix m_Y; // curves.train in R
        std::vector<int> m_reconst_fcts; //vector of indices of curves to reconstruct 
        double m_alpha;
        int m_K; //check type
        NumericVector m_t_points; //check type
        int m_maxBins;
        int m_nRegGrid;
};


class ReconstructionKraus : public ReconstructionBase{
    public:
        //constructor for Kraus
        ReconstructionKraus(const NumericMatrix& Y, double alpha, int K, 
                            NumericVector& t_points, int nRegGrid, int maxBins) : 
                            ReconstructionBase(Y, alpha, K, t_points, nRegGrid, maxBins) //per unique_ptr
        void reconstructCurves() override;//metodo che sarà chiamato da R

        const std::vector<double>& meanKraus(const NumericMatrix&);//oppure return void?
        const std::vector<std::vector<double>>& covKraus(const NumericMatrix&); //oppure return void?
        

        // getter e setter mean e cov
        NumericVector mean() const { return wrap(m_mean); } //returns value computed by meanKraus
        NumericMatrix cov() const; //returns value computed by covKraus

        double gcvKraus(const std::vector<std::vector<double>>& covMat, const std::vector<double>& meanVec, 
                        const NumericMatrix& X, const bool M_bool_vec, const Numeric& alpha) const;
        TUPLA reconstKraus_fun(const std::vector<std::vector<double>>& covMat, 
                               const std::vector<double>& X_cent_vec, 
                               const Numeric& alpha); //ritorna una lista => how to deal?

        //override
        double alpha() const override {return m_alpha; }
        int K() const override { Rcout << "no K for Kraus method" << std::endl; }
        NumericVector t_points() const override { Rcout << "no t_points for Kraus method" << std::endl; }
        int max_bins() const override { Rcout << "no bins for Kraus" << std::endl; }
        int n_reg_grid() const override { Rcout << "no n_reg_grid for Kraus" << std::endl; }

    private:
        std::vector<double> m_mean; // set to NA default? capire se mettere return type NumericVector
        std::vector<std::vector<double>> m_cov; //default? return type NumericMatrix?

};

class ReconstructionKLAl : public ReconstructionBase{//capisci se K era un double o un int in R
    public: 
        ReconstructionKLAl(const NumericMatrix& Y, double alpha = 0.0, int K = 0, 
                           NumericVector& t_points = NumericVector(), int nRegGrid = 0, int maxBins = 0) : 
                           ReconstructionBase(Y, alpha, K, t_points, nRegGrid, maxBins) //per unique_ptr
        void reconstructCurves() override;

        int K() const override { return m_K; }
        NumericVector t_points() const override { return m_t_points; }
        int max_bins() const override { return m_maxBins; }
        int n_reg_grid() const override { return m_nRegGrid; }
        double alpha() const override { Rcout<< "no alpha for KL method" << std::endl; }

    private:

};

#endif RECONSTRUCTION_H