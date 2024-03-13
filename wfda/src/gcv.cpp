#include "gcv.h"
//correggi la funzione -> X_gcv e accesso a m_YY
double gcv::value(double x){//gcvKrausFun
    int n = m_YY.ncol();
    NumericVector rss = NumericVector::create(n);
    double df = 0.0;
    for(int j=0; j < n ; ++j){
        NumericVector X_gcv = m_YY(_,j);//c'è qualcosa di sbagliato perchè non è chiamato da reconstKraus_fun
        //potrei modificare man mano il data member m_YY
        X_gcv[m_M_bool] = NA_REAL; //secondo me non posso accedere così, cambia! o forse si
        List result_tmp = reconstKraus_fun(m_YY, m_mean_gcv, m_cov_gcv, m_index, x);//x = alpha
        //devo cambiare qualcosa in Kraus perchè così non gli arriva X_gcv mascherato
        //magari farlo per il caso in cui è chiamato da gcvKraus e quando no
        //fare un overload??? per esempio gli passo un termine in più ed è bool
        NumericVector X_fit = Numericvector(result_tmp["X_cent_reconst_vec"].begin(),result_tmp["X_cent_reconst_vec"].end()) + m_mean_gcv;//if m_mean is vector<double> need to use as<NumericVector>(m_mean)
        rss[j] = sum(pow((X_fit[m_M_bool] - m_YY[m_M_bool,j]),2));//è sbagliato l'accesso alla matrice, posso usare LogicalVector solo su Vector

        if(j == n-1)
            df = result_tmp["df"];
    }

    double gcv = sum(rss)/(std::pow(1-df/n,2));
    return gcv;
}
