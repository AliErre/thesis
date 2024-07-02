#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H
#include <vector>
#include <list>
#include <RcppArmadillo.h>

#include <omp.h>
using namespace Rcpp;

class Genetic{
    public:
    //fare anche altri costruttori più semplici per overload
        //aggiornare i costruttori perchè manca l'inizializzazione di alcuni data member
        Genetic(const List& xlist, const List& blistdefault, const arma::mat& curves, const List& curvesfd, const arma::vec& tpoints, const arma::uvec& events):
                m_xlist(xlist), m_b_list_def(blistdefault), m_curves(curves), m_curvesfd(curvesfd), m_tpoints(tpoints), m_events(events) {m_n_par = xlist.length();}
        void multistart();
        std::pair<arma::vec,arma::vec> pwmse(List&);//move from member to free function, but add arguments

        arma::mat get_P(int );
        arma::vec get_V(int );
        arma::mat get_Pin(int );
        arma::vec get_Vin(int );
        List get_blist(int );

    private:
        Environment base_env = Environment::base_env();
        Function set_seed = base_env["set.seed"];
        size_t m_niter = 5;  //era 15
        size_t m_B = 10;  //era 10
        size_t m_P_size = 20; //era 20
        
        double m_r = 0.5;
        arma::vec perturbation_vec{0.5, 2};
        size_t m_n_multistart = 5;
        int m_seed_in = 22;
        List m_xlist;
        size_t m_n_par= m_xlist.length();
        List m_b_list_def;
        arma::mat m_curves;
        List m_curvesfd;
        arma::vec m_tpoints;
        arma::uvec m_events;
        Environment fda = Environment::namespace_env("fda");
        Function create_constant_basis = fda["create.constant.basis"];
        Function fd = fda["fd"];
        Function eval_fd = fda["eval.fd"];
        Function fdPar = fda["fdPar"];
        Function inprod = fda["inprod"];
        Function eval_penalty = fda["eval.penalty"];
        Function smooth_basis = fda["smooth.basis"];
        Function int2Lfd = fda["int2Lfd"];
        Function times_fd = fda["times.fd"];
        Function sum_fd = fda["sum.fd"];

        std::map<int, arma::mat> m_P_initial;//save(P.intial, V.initial)
        std::map<int, arma::vec> m_v_initial;//save(P.intial, V.initial)
        std::map<int, arma::mat> m_P;//save(P, V)
        std::map<int, arma::vec> m_v;//save(P, V)
        //fare getter del tipo get_data(int indice){check esistenza indice}
        std::map<int, List> m_blist_eaas;//save(blist.eaas)
        //dopo i test fare il getter ed esporli con .property()
};
#endif //GENETIC_ALGORITHM_H