#ifndef GENCROSSVAL_H
#define GENCROSSVAL_H
#include "optimize.h"
#include "reconstKraus.h"

class gcv : public Optim{
    public:
        gcv(NumericMatrix& Y, NumericVector& mean_vec,
            NumericMatrix& cov_mat,unsigned index) :
            m_YY(Y), m_cov_gcv(cov_mat), m_mean_gcv(mean_vec), m_index(index) {};
        double value(double x) override;
    private:
        NumericMatrix m_YY;//X_compl_mat
        NumericVector m_mean_gcv;//m_mean
        NumericMatrix m_cov_gcv;//m_cov
        LogicalVector m_M_bool;
        unsigned m_index;       

};


#endif GENCROSSVAL_H