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
        ReconstructionBase(const std::string& method, const NumericMatrix& Y) 
                          : m_method(method), m_Y(Y) 
                          { m_reconst_fcts = find_obs_inc(Y); meanRows(); covMatrix();}

        // specialize
        virtual ~ReconstructionBase() = default; // polymorphism => need virtual destructor
        virtual List reconstructCurve(Nullable<double>, bool, const Nullable<NumericVector>&, Nullable<int>, Nullable<int>, Nullable<int>) = 0; //i,alpha,K,t_points,nRegGrid,maxBins 
        //devo aggiungere i default
        // same for all derived
        std::vector<int> find_obs_inc(const NumericMatrix&) const; //farne una free function?
        IntegerVector reconst_fcts() const;//getter
        const NumericVector& meanRows(); //farne una free function?mettere return type NumericVector? penso sia meglio!!!
        const NumericMatrix& covMatrix(); //farne una free function?    
        const NumericVector& mean() const { return m_mean; } //beware wrap() copies the object -> heavy when data is big
        const NumericMatrix& cov() const {return m_cov;} //returns value computed by covKraus
        //since it's needed both for klal and klnoal and I dont want to duplicate code
        

    protected:
        std::string m_method;
        NumericMatrix m_Y; // curves.train in R
        std::vector<int> m_reconst_fcts; //vector of indices of curves to reconstruct 
        NumericVector m_mean; //capire se mettere return type NumericVector
        NumericMatrix m_cov; //return type NumericMatrix cause it can contain NA
};


class ReconstructionKraus : public ReconstructionBase{
    public:
        //constructor for Kraus
        ReconstructionKraus(const std::string& method, const NumericMatrix& Y) : ReconstructionBase(method, Y) {};

        //override
        List reconstructCurve(Nullable<double> alpha, bool all, const Nullable<NumericVector>& periods_nullable, Nullable<int> K, Nullable<int> maxBins, Nullable<int> nRegGrid_nullable) override;//metodo che sarà chiamato da R

};

class ReconstructionExtrapolation : public ReconstructionBase{
    public:
        //constructor for Kraus
        ReconstructionExtrapolation(const std::string& method, const NumericMatrix& Y) : ReconstructionBase(method, Y) {};

        //override
        List reconstructCurve(Nullable<double>, bool, const Nullable<NumericVector>&, Nullable<int>, Nullable<int>, Nullable<int>) override;//metodo che sarà chiamato da R
        //NumericVector è T.periods in R

};

class ReconstructionKL : public ReconstructionBase{//capisci se K era un double o un int in R
    public: 
        ReconstructionKL(const std::string& method, const NumericMatrix& Y) : ReconstructionBase(method, Y){}; 
        List reconstructCurve(Nullable<double>, bool, const Nullable<NumericVector>&, Nullable<int>, Nullable<int>, Nullable<int>) override;
        //NumericVector è t.points, K, nRegGrid, maxBins
        //Ly, Lu, reconst_fcts,CEscores, center, maxBins 
        void myfpca(std::vector<std::vector<double>>&, const std::vector<std::vector<double>>&, 
                    bool, bool, Nullable<int>, bool);
    private: 
        std::pair<std::vector<double>, NumericMatrix> m_Y_preprocessed;
        std::vector<std::vector<double>> m_observed_period;
        arma::mat m_cov_est;
        NumericVector m_mu;
        std::vector<arma::vec> m_muO;
        std::vector<double> m_scoresO;
        List m_CE_scoresO;
        arma::mat m_efunctions;
        std::vector<arma::mat> m_efunctionsO;
        std::vector<NumericMatrix> m_efun_reconst;
        arma::vec m_evalues;
        std::vector<arma::vec> m_evaluesOO;
        std::vector<arma::vec> m_obs_argvalsO;
        std::vector<arma::uvec> m_locO;
        double m_sigma2;        

};


#endif