#ifndef GENCROSSVAL_H
#define GENCROSSVAL_H
#include "optimize.h"
#include "reconstKraus.h"

class gcv : public Optim{
    public:
        gcv(const NumericMatrix Y,const NumericVector& mean_vec,
            const NumericMatrix& cov_mat) :
            m_YY(Y),  m_mean_gcv(mean_vec), m_cov_gcv(cov_mat){}
        double value(double x) override;
        void set_bool(const LogicalVector& m){ m_M_bool = m;}
    private:
        NumericMatrix m_YY;//X_compl_mat
        NumericVector m_mean_gcv;//m_mean
        NumericMatrix m_cov_gcv;//m_cov
        LogicalVector m_M_bool;
        int m_index;//mi sa che non è usato, togli       

};


#endif