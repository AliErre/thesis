#include "gcv.h"
//correggi la funzione -> aggiungi metodi per l'update di m_YY, così non ricreo l'oggetto gcv ogni volta 
double gcv::value(double x){//gcvKrausFun
    int n = m_YY.ncol();
    NumericVector rss = NumericVector::create(n);
    double df = 0.0;
    for(int j=0; j < n ; ++j){
        //maschera valori in m_YY(_,j)
        NumericVector original(sum(m_M_bool));
        for(int r = 0;r < m_M_bool.length() ;++r){
            if(m_M_bool[r]) 
            {
                original.push_back(m_YY(r,j));
                m_YY(r,j) = NA_REAL;//nel codice è il ruolo di X_gcv
            } 
        }
        List result = reconstKraus_fun(m_YY, m_mean_gcv, m_cov_gcv, j, x);//x = alpha
        NumericVector X_fit = result["X_cent_reconst_vec"] + m_mean_gcv;
        NumericVector X_fit_subset = X_fit[m_M_bool];//non potevo fare la chiamata direttamente con X_fit[m_M_bool]
        rss[j] = sum(pow(X_fit_subset-original,2));//forse conviene scrivere la nostra pow(,), perchè se non ricordo male pow() non è efficiente

        if(j == n-1)
            df = result["df"];
    }

    double gcv = sum(rss)/(std::pow(1-df/n,2)); //in realtà questo std::pow() non è efficiente, me lo devo fare io
    return gcv;
}
