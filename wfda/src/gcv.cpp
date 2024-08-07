#include "gcv.h" 
double gcv::value(double x){//gcvKrausFun
    int n = m_YY.ncol();//X_Compl_mat
    int bl = m_M_bool.length();
    NumericVector rss = NumericVector(n);
    double df = 0.0;
    for(int j=0; j < n ; ++j){
        //maschera valori in m_YY(_,j)
        NumericVector original(sum(m_M_bool));
        int index = 0;
        for(int r = 0;r < bl ;++r){
            if(m_M_bool[r]) 
            {
                original[index++] = m_YY(r,j);
                m_YY(r,j) = NA_REAL;//nel codice è il ruolo di  X_gcv[M_bool_vec] <- NA
            } 
        }
        std::tuple<NumericVector, double, double, arma::vec, arma::uvec> result = reconstKraus_fun(m_YY, m_mean_gcv, m_cov_gcv, j, x);//x = alpha
        //reset
        index = 0;
        for(int r = 0;r < bl ;++r){
            if(m_M_bool[r]) 
            {
                m_YY(r,j) = original[index++];
            } 
        }
        NumericVector X_cent_reconst_vec = std::get<0>(result);
        NumericVector X_fit = X_cent_reconst_vec + m_mean_gcv;
        NumericVector X_fit_subset = X_fit[m_M_bool];//non potevo fare la chiamata direttamente con X_fit[m_M_bool]
        for(int i = 0;i < original.length();++i)
            rss[j] = sum(pow(X_fit_subset-original,2));//pow() non è efficiente

        if(j == n-1)
            {df = std::get<1>(result);}
    }
    double gcv = sum(rss)/((1 - df / n) * (1 - df / n)); // std::pow() non è efficiente, farlo io
   return gcv;
}