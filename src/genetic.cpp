#include "genetic.h"
#include "helper_genetic.h"
#include "sample.h"

std::pair<arma::vec, arma::vec> Genetic::pwmse(List& blist)//toglierla dai member e aggiungere parametri
{   //non credo di volerla togliere dai method perchè così può prendersi le funzioni messe come data member
   // size_t n = m_curves.n_cols;
    size_t N = m_curves.n_rows;
    size_t p = m_xlist.length();
    //training and test set split
    set_seed(140996);
    arma::uvec unique_events = arma::unique(m_events);
    size_t n_events = unique_events.size();
    size_t n_test_events = static_cast<size_t>(n_events / m_B);
    arma::uvec shuffled_unique_events = sample(unique_events, n_events, false);
    arma::mat mse_eval(N, m_B, arma::fill::zeros);
    arma::mat sigma_eval(N, m_B, arma::fill::zeros);
    //arma::mat mse_ita18
    //arma::mat sigma_ita18
    for(size_t b = 1; b <= m_B; ++b)
    {
        arma::uvec idxs = arma::conv_to<arma::uvec>::from(arma::regspace((b - 1)*n_test_events +1, std::min(b * n_test_events, n_events)));
        //indice di posizione con cui accedere a shuffled_unique_events
        arma::uvec test_events = shuffled_unique_events(idxs-1);
        arma::uvec test_indices;
       for(const auto& t:test_events)
        {
            arma::uvec indices = arma::find(m_events == t);
            test_indices = arma::join_cols(test_indices, indices);
        }
        arma::mat curves_test = m_curves.cols(test_indices);
        arma::uvec train_indices = arma::conv_to<arma::uvec>::from(arma::regspace(1, m_events.size())-1);
        train_indices.shed_rows(test_indices);
        
        arma::mat curves_train(m_curves.n_rows, train_indices.size());
        curves_train = m_curves.cols(train_indices);
         
        /*arma::uvec train_indices(train_indices_.size());
        size_t column = 0;
        for(const auto& s:train_indices_)
        {
            train_indices[column++] = s;
        }*/
        List curves_fd_train = clone(m_curvesfd);
        arma::mat curves_fd_coefs = curves_fd_train["coefs"];
        arma::mat curves_fd_train_coefs = curves_fd_coefs.cols(train_indices);
        curves_fd_train["coefs"]= curves_fd_train_coefs;
        /*wgts.fd.train <- wgts.fd # no questo
        wgts.fd.train$coefs <- wgts.fd$coefs[,train]
    
        data.test      <- data.f[test,]# default NULL. Salta
        data.train     <- data.f[train,]*/
        List xlist_test(m_xlist.length());
        List xlist_train(m_xlist.length());
        for(size_t i = 0; i < p; i++)
        {
            if (Rf_inherits(m_xlist[i], "fd")) {
                
                List xlist_i = m_xlist[i];
                arma::mat coefs = xlist_i["coefs"];//vedere se il casting da NumericMatrix ad arma::mat è automatico
                arma::mat coefs_test = coefs.cols(test_indices);
                arma::mat coefs_train = coefs.cols(train_indices);

                List x_test_i = clone(xlist_i);
                x_test_i["coefs"] = coefs_test;
                xlist_test[i] = x_test_i;

                List x_train_i = clone(xlist_i);
                x_train_i["coefs"] = coefs_train;
                xlist_train[i] = x_train_i;
            }else if(Rf_isNumeric(m_xlist[i])) {
                
                arma::vec xlist_i = as<arma::vec>(m_xlist[i]);
                arma::vec xlist_i_test = xlist_i(test_indices);
                arma::vec xlist_i_train = xlist_i(train_indices);
                NumericVector xlist_i_test_n = wrap(xlist_i_test);
                NumericVector xlist_i_train_n = wrap(xlist_i_train);//perchè dopo my_fRegress_argcheck ha i controlli come: Rf_isNumeric
                xlist_test[i] = xlist_i_test_n;//se poi così è giusto, fai le prove mettendo direttamente xlist_i_test
                xlist_train[i] = xlist_i_train_n;//prova mettendo direttamente xlist_i_train
            }else if(Rf_isMatrix(m_xlist[i])){
                 
                arma::mat xlist_i = as<arma::mat>(m_xlist[i]);
                arma::mat rows = xlist_i.rows(test_indices);
                arma::vec row_vec = rows.col(0);
                NumericVector row = wrap(row_vec);
                xlist_test[i] = row;//fai le prove mettendo row_vec
                rows = xlist_i.rows(train_indices);
                row_vec = rows.col(0);
                row = wrap(row_vec);
                xlist_test[i] = row;  //fai le prove mettendo row_vec              
            }
        }
        std::tuple<List, List, List, List, List, arma::mat, arma::vec> mod;
        //qua c'era if else wgtsaflag
        mod =  weighted_fRegress(curves_fd_train, xlist_train, blist, R_NilValue, false, create_constant_basis, fd, 
                                 inprod, eval_penalty, eval_fd, smooth_basis, times_fd, sum_fd, fdPar, int2Lfd);
       
        NumericVector range= {arma::min(m_tpoints), arma::max(m_tpoints)};
        List onebasis = create_constant_basis(range);//fai un test, magari è da fare un wrap
        List onesfd = fd(1, onebasis);
        List xfdlist_test = xlist_test;
        for(size_t i = 0; i < p; i++)
        {
            if(Rf_isNumeric(xfdlist_test[i]))
            {
                arma::vec xfdlist_i = xfdlist_test[i];
                arma::mat reshaped = arma::reshape(xfdlist_i, 1, test_indices.size());
                xfdlist_test[i] = fd(reshaped, onebasis);
            } 
        }
        List curves_hat = predict_fRegress(mod, xfdlist_test, m_tpoints, eval_fd, smooth_basis);
        arma::mat curves_hat_v = as<arma::mat>(eval_fd(_["evalarg"] = m_tpoints, _["fdobj"] = curves_hat));//check se puoi
        arma::vec nT(N);
        for(size_t t = 0; t < N; t++)
        {
            if(curves_test.n_rows > 1)
            {
                arma::uvec no_na = arma::find_finite(curves_test.row(t));
                nT[t] = no_na.size();
            }else{
                nT[t] = 1;
            }
        }
        arma::mat diff = curves_test - curves_hat_v;
        arma::mat diff2 = arma::square(diff);
        arma::uvec is_na = arma::find_nonfinite(diff2);
        diff2.elem(is_na).zeros();
        arma::vec rowsums = arma::sum(diff2, 1);
        mse_eval.col(b-1) = rowsums/nT;
        sigma_eval.col(b-1) = arma::pow(rowsums/(nT - p), 0.5);
    }

    arma::vec mse_t = arma::sum(mse_eval, 1)/m_B;
    arma::vec sigma_t = arma::sum(sigma_eval, 1)/m_B;
    return std::make_pair(mse_t, sigma_t);
}

void Genetic::multistart(){//finirla: mancano tutte le cose che salvava nei file, da settare come data member
    arma::mat Pinitial(m_P_size, m_n_par);
  
    //for(size_t b = 1; b <= m_n_multistart; b++ ){
    int b=1;
    int seed_b = m_seed_in + b;  // per avere seed come R
        set_seed(seed_b);

        arma::vec exp_range = arma::regspace(-5,2);
        arma::vec tens(exp_range.size(), arma::fill::value(10));
        arma::vec lambda_range = arma::pow(tens, exp_range);
        for(size_t i = 0; i < Pinitial.n_rows; i++){
            arma::vec result = sample(lambda_range, m_n_par, true);
            Pinitial.row(i) = result.t();
        }
        arma::vec v_initial(m_P_size);
        List blist;
        for(size_t i = 0; i < m_P_size; i++){
            blist = clone(m_b_list_def);
            for(size_t reg = 0; reg < m_n_par; reg++)
            {
                List b = blist[reg];
                b["lambda"] = Pinitial(i, reg);
                blist[reg] = b;
            }
            std::pair<arma::vec, arma::vec> ret = pwmse(blist);//b-fold
            arma::vec mse = ret.first;
            double sum = arma::accu(arma::square(mse));
            v_initial[i] =sum/static_cast<double>(mse.size());
        }//save(P.initial, V.initial, file=paste0("data/calibration/Initial_blist_select_seed",seed,".RData"))
        m_P_initial[b + 1] = Pinitial;//to have R indexing
        m_v_initial[b + 1] = v_initial;
        arma::vec v = v_initial;
        arma::mat P = Pinitial;
        arma::vec v_sort;
        arma::uvec indx_in;
        size_t size = static_cast<size_t>(m_r * m_P_size);
        arma::mat Q;
        arma::mat Q_prime;
        arma::vec z;
        arma::vec perturbation;
        arma::mat pert_mat(m_n_par, size);//because filled column major
        arma::vec z_prime(size);
        for(size_t k = 0; k < m_niter; k++)
        {
            v_sort = arma::sort(v, "descend");
            arma::vec v_sort_sub = v_sort.subvec(0, size-1);
            arma::vec a(v.size(), arma::fill::zeros);
            for (auto i:v_sort_sub) a = a +  (v==i);//which(v %in% v_sort_subvec[size])
            indx_in = arma::find(a);//finds non zero elements. non zero if it occured in subvec
            if(indx_in.size() > size){//avoid ties
                set_seed(seed_b); 
                indx_in = sample(indx_in, size, false);
            }
            Rcout<<"ind "<<k;
            //Q.set_size(indx_in.size(), P.n_cols);
            
            Q = P.rows(indx_in);//forse non serve creare P, va bene P_initial
            z.set_size(indx_in.size());
            z = v(indx_in);
            perturbation.set_size(size * m_n_par);
            set_seed(seed_b);  // per avere garanzia sia lo stesso
            perturbation = sample(perturbation_vec, size * m_n_par, true);
            pert_mat = arma::reshape(perturbation, size, m_n_par);
            Q_prime = Q % pert_mat;   // elem per elem
            arma::uvec indx_ext = arma::find(Q_prime < 1e-5 || Q_prime > 1e2);//matrices are column major in arma
            Q_prime(indx_ext) = sample(lambda_range, indx_ext.size(), true);

            for(size_t i = 0; i < size; i++)
            {
                blist = clone(m_b_list_def);
                for(size_t reg = 0; reg < m_n_par; reg++)
                {
                    List b = blist[reg];
                    b["lambda"] = Q_prime(i, reg);
                    blist[reg] = b;
                }            
                std::pair<arma::vec, arma::vec> ret = pwmse(blist);//b-fold se la lascio data member e m_B è lo stesso in entrambe le chiamate posso togliere m_B dagli argomenti           
                arma::vec mse = ret.first;
                double sum = arma::accu(arma::square(mse)); 
                z_prime[i] = sum/mse.size();  
            }
            arma::mat P_in = P.rows(indx_in);
            P = arma::join_cols(P_in, Q_prime);
            arma::vec v_in = v(indx_in);
            v = arma::join_vert(v_in, z_prime);
        }//save(P, V, file=paste0("data/calibration/P_V_blist_select_seed",seed,".RData"))
        
        m_P[b + 1] = P;
        m_v[b + 1] = v;
        arma::uword best = v.index_min();//da per scontato sia uno
        arma::rowvec lambda_f = P.row(best);
        arma::vec lambda_opt(P.n_cols);
        lambda_opt = lambda_f.t();
        for(size_t reg = 0; reg < m_n_par; ++reg)
        {
            List b = m_b_list_def[reg];
            b["lambda"] = lambda_opt[reg];
            m_b_list_def[reg] = b;
        }//save(blist, file=paste0('data/calibration/blist_EAASS_seed',seed,'.RData'))
        m_blist_eaas[b + 1] = m_b_list_def;
      
    }  // for } 


arma::mat Genetic::get_P(int elem){
    return m_P[elem];
}
arma::vec Genetic::get_V(int elem){
    return m_v[elem];
}
arma::mat Genetic::get_Pin(int elem){
    return m_P_initial[elem];
}
arma::vec Genetic::get_Vin(int elem){
    return m_v_initial[elem];
}
List Genetic::get_blist(int elem){
    return m_blist_eaas[elem];
}