#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H
#include <memory>
#include <vector>
#include "optimize.h"
#include "reconstKraus.h"
using namespace Rcpp;

class ReconstructionBase {
    public:
    //nota: references point the memory area of the R object!!
        ReconstructionBase(const NumericMatrix& Y) 
                          : m_Y(Y)
                          { m_reconst_fcts = find_obs_inc(Y); meanRows(); covMatrix();}

        // specialize
        virtual ~ReconstructionBase() = default; // polymorphism => need virtual destructor
        virtual List reconstructCurve(double, bool, const NumericVector&, int, int) = 0; //i,alpha,K,t_points,nRegGrid,maxBins 

        // same for all derived
        std::vector<int> find_obs_inc(const NumericMatrix&) const; //farne una free function?
        IntegerVector reconst_fcts() const;//getter
        const NumericVector& meanRows(); //farne una free function?mettere return type NumericVector? penso sia meglio!!!
        const NumericMatrix& covMatrix(); //farne una free function?    
        // getter e setter mean e cov
        NumericVector mean() const { return m_mean; } //beware wrap() copies the object -> heavy when data is big
        NumericMatrix cov() const {return m_cov;} //returns value computed by covKraus

    protected:
        NumericMatrix m_Y; // curves.train in R
        std::vector<int> m_reconst_fcts; //vector of indices of curves to reconstruct 
        NumericVector m_mean; //capire se mettere return type NumericVector
        NumericMatrix m_cov; //return type NumericMatrix cause it can contain NA
};


class ReconstructionKraus : public ReconstructionBase{
    public:
        //constructor for Kraus
        ReconstructionKraus(const NumericMatrix& Y) : ReconstructionBase(Y) {};

        //override
        List reconstructCurve(double, bool, const NumericVector&, int, int) override;//metodo che sarà chiamato da R

};

class ReconstructionExtrapolation : public ReconstructionBase{
    public:
        //constructor for Kraus
        ReconstructionExtrapolation(const NumericMatrix& Y) : ReconstructionBase(Y) {};

        //override
        List reconstructCurve(double, bool, const NumericVector&, int, int) override;//metodo che sarà chiamato da R
        //NumericVector è T.periods in R

};

class ReconstructionKLAl : public ReconstructionBase{//capisci se K era un double o un int in R
    public: 
        ReconstructionKLAl(const NumericMatrix& Y): //dovrò aggiungere questi argomenti a reco_factory
                           ReconstructionBase(Y) {}; 
        List reconstructCurve(double, bool, const NumericVector&, int, int) override;
        //NumericVector è t.points, K, nRegGrid, maxBins
        void myfpca(const std::vector<std::vector<double>>&, const std::vector<std::vector<double>>&, bool, bool, int, bool);
        //Ly, Lu, reconst_fcts,CEscores, center, maxBins 

    private:
        std::pair<std::vector<double>, NumericMatrix> m_Y_preprocessed;
        std::vector<std::vector<double>> m_observed_period;
        arma::mat m_cov_est;
        NumericVector m_mu;
        List m_muO;
        List m_scoresO;
        List m_CE_scoresO;
        arma::mat m_efunctions;
        List m_efunctionsO;
        List m_efun_reconst;
        arma::vec m_evalues;
        List m_evaluesOO;
        List m_obs_argvalsO;
        List m_locO;
        double m_sigma2;
};

#endif