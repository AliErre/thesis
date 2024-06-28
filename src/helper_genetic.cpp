#include "helper_genetic.h"
bool is_fdPar(const List& obj) {
  if (obj.containsElementNamed("fd")) {
    List fdobj = obj["fd"];
    return fdobj.containsElementNamed("basis") && fdobj.containsElementNamed("coefs");  // all elements needed to be an fdPar
  }
  return false;
}


void my_fRegressArgCheck(List& yfd, List& xfdlist, List& betalist, const Function& create_constant_basis, const Function& fd, const Function& fdPar)
{
     //Environment fda = Environment::namespace_env("fda");
     //Function fdPar = fda["fdPar"];//aggiungere alle function passate, metti la reference
    if(!(is_fdPar(yfd) || Rf_inherits(yfd, "fd") || Rf_isNumeric(yfd) || Rf_isMatrix(yfd))){
        stop("First argument is not of class 'fdPar', 'fd', 'numeric' or 'matrix'.");
    }
    if(is_fdPar(yfd)){ yfd = yfd["fd"];}//siccome passo le reference mi sa che posso fare void e non ritornare arglist in pwmse

    int N = 0;
    NumericVector rangeval;
    if(Rf_inherits(yfd, "fd")){ 
        NumericMatrix coefs = yfd["coefs"]; 
        N = coefs.ncol();
        List yfdbasis = yfd["basis"];
        rangeval = yfdbasis["rangeval"];
    }else{
        N = Rcpp::as<Rcpp::NumericVector>(yfd).size();
        rangeval = {0, 1};
    }
    
    int p = xfdlist.length();
    if(betalist.length() != p){ stop("number of regression coefficients does not match number of independent variables");}
    
    //xfdlist: if the object is a vector of length N, it is converted to a functional data object with a constant basis
    List onebasis = create_constant_basis(_["rangeval"] = rangeval);
    List onesfd = fd(1, onebasis);
    for(int j = 0; j < p; j++)
    {
        if(Rf_isNumeric(xfdlist[j])){
            NumericVector xfdj = as<NumericVector>(xfdlist[j]);
            NumericMatrix mat(1, N, xfdj.begin());
            xfdlist[j] = fd(mat, onebasis);
        }
        if(!(Rf_inherits(xfdlist[j],"fd") || Rf_isNumeric(xfdlist[j]) || Rf_isMatrix(xfdlist[j]))){stop("error in xfdlist");}
    }

    for(int j = 0; j < p; j++)
    {
        List betafdParj = betalist[j];
        if(Rf_inherits(betafdParj, "fd") || Rf_inherits(betafdParj, "basisfd")){
            betafdParj = fdPar(betafdParj);
            betalist[j] = betafdParj;
        }
        if(!Rf_inherits(betafdParj, "fdPar")){
            stop("betalistj is not a fdpar object");
        }
    }
}


List predict_fRegress(const std::tuple<List, List, List, List, List, arma::mat, arma::vec>& mod, const List& xlist, const arma::vec& tpoints,
                      const Function& eval_fd, const Function& smooth_basis)
{
    /*Environment fda = Environment::namespace_env("fda");
    Function eval_fd =  fda["eval.fd"];
    Function smooth_basis = fda["smooth.basis"];*/
    List y = std::get<0>(mod);
    List betaestlist = std::get<3>(mod);
    size_t p = static_cast<size_t>(xlist.length());
    List ybasisobj = y["basis"];
    List xlist_first = xlist[0];
    arma::mat coefs = as<arma::mat>(xlist_first["coefs"]);
    int n = coefs.n_cols;
    int ynbasis = ybasisobj["nbasis"];
    int nfine = std::max(501, 10 * ynbasis + 1);
    arma::vec tfine = arma::linspace(arma::min(tpoints), arma::max(tpoints), nfine);
    arma::mat yhatmat(nfine, n, arma::fill::zeros);
    for(size_t j = 0; j < p; j++)
    {
        arma::mat xfdj = as<arma::mat>(eval_fd(_["evalarg"] = tfine, _["fdobj"] = xlist[j]));//controlla in R se veramente da una matrice
        List betafdParj = betaestlist[j];
        List betafdj = betafdParj["fd"];
        arma::mat betavecj = as<arma::mat>(eval_fd(_["evalarg"] = tfine, _["fdobj"] = betafdj));
        yhatmat = yhatmat + xfdj.each_col() % arma::vectorise(betavecj);
    }
    List res = smooth_basis(_["argvals"] = tfine, _["y"] = yhatmat, _["fdParobj"] = ybasisobj);
    List yhatfdobj = res["fd"];
    return yhatfdobj;
}



std::tuple<List, List, List, List, List, arma::mat, arma::vec> 
weighted_fRegress(List& y, List& xfdlist, List& betalist, const Nullable<List>& wgts, bool weighted,
                  const Function& create_constant_basis, const Function& fd, const Function& inprod,
                  const Function& eval_penalty, const Function& eval_fd, const Function& smooth_basis,
                  const Function& times_fd, const Function& sum, const Function& fdPar, const Function& int2Lfd)
{
    //vedi se ha senso mettere le reference
    //così vengono cambiate nella funzione
    /*Environment fda("package:fda");        // import all fda functions used in this piece of code
    Function create_constant_basis = fda["create.constant.basis"];
    Function fd = fda["fd"];
    Function inprod = fda["inprod"];
    Function eval_penalty = fda["eval.penalty"];
    Function eval_fd = fda["eval.fd"];
    Function int2Lfd = fda["int2Lfd"];
    Function smooth_basis = fda["smooth.basis"];
    Function times_fd = fda["times.fd"];
    Function sum = fda["sum.fd"];
    Function fdPar = fda["fdPar"];*/
    if(is_fdPar(y)){
     List yfdobj = y["fd"];
     y = yfdobj;
    }

    my_fRegressArgCheck(y, xfdlist, betalist, create_constant_basis, fd, fdPar);
    List yfdobj = y;//["yfd"];
    //xfdlist = std::get<1>(arglist);//["xfdlist"];
    //betalist = std::get<2>(arglist);//["betalist"];  // update my data, xfdlist deve essere riempito con dati fd  
    size_t p = static_cast<size_t>(xfdlist.length());   
    arma::mat ycoef = yfdobj["coefs"];
    int N = ycoef.n_cols;
    List ybasisobj = yfdobj["basis"];
    NumericVector range = ybasisobj["rangeval"];
    int nbasis = ybasisobj["nbasis"];
    List onesbasis = create_constant_basis(range);//vedi se serve un wrap()
    List onesfd= fd(1, onesbasis); 

    int ncoef = 0;     //calcolo numero coeff in regression
    for (size_t j = 0; j < p; j++) {
        List betafdParj = betalist[j];
        if (betafdParj["estimate"]) {
            List fdobj = betafdParj["fd"];
            List basisobj = fdobj["basis"];
            int ncoefj = basisobj["nbasis"];   // senza spacchettare così non funziona
            ncoef += ncoefj;
        }
    }
 
    arma::mat Cmat(ncoef, ncoef,arma::fill::zeros);
    arma::vec Dmat(ncoef, arma::fill::zeros);
    int mj2 = 0;
    int mj1 = 0;
    for(size_t j = 0; j < p; j++)
    {
        List betafdParj = betalist[j];
        if(betafdParj["estimate"])
        {
            List betafdj = betafdParj["fd"];
            List betabasisj = betafdj["basis"];
            int ncoefj = betabasisj["nbasis"];
            //row indices
            mj1 = mj2 + 1;
            mj2 += ncoefj;
            arma::uvec indexj = arma::conv_to<arma::uvec>::from(arma::regspace(mj1, mj2));
            //compute right side of equation DMAT
            //compute weight function for DMAT
            List xfdj = xfdlist[j];
            List xyfdj;
            List wgts_;
            if(wgts.isNull())
            {
                xyfdj = times_fd(xfdj, yfdobj);
            }else{
                wgts_ = wgts.get();
                List w = times_fd(xfdj, wgts_);
                xyfdj = times_fd(w, yfdobj);
            }
            List wtfdj = sum(xyfdj, fd);
            
            //jth component of DMAT
             NumericMatrix Cmatjk = inprod(_["fdobj1"] = betabasisj,_["fdobj2"] = onesfd, _["Lfdobj1"] = int2Lfd(0), _["Lfdobj2"] = int2Lfd(0), _["rng"] = range, _["wtfd"] = wtfdj);
             arma::mat Dmatj = as<arma::mat>(Cmatjk);
            Dmat(indexj-1) = Dmatj; 
            int mk1 = 0;
            int mk2 = 0;
            for(size_t k = 0; k <= j; k++)
            {
                List betafdPark = betalist[k];
                if(betafdPark["estimate"]){
                    List betafdk = betafdPark["fd"];
                    List betabasisk = betafdk["basis"];
                    int ncoefk = betabasisk["nbasis"];
                    mk1 = mk2 + 1;
                    mk2 += ncoefk;
                    arma::uvec indexk = arma::conv_to<arma::uvec>::from(arma::regspace(mk1, mk2));
                    //set up the weight function for CMAT
                    List xfdk = xfdlist[k];
                    List xxfdjk;
                    if(wgts.isNull())
                    {
                        xxfdjk = times_fd(xfdj, xfdk);
                    }else{
                        List ww = times_fd(xfdj, wgts_);
                        List xxfdjk = times_fd(ww, xfdk);
                    }
                    List wtfdjk = sum(xxfdjk, fd);
                    arma::mat Cmatjk = as<arma::mat>(inprod(_["fdobj1"] = betabasisj, _["fdobj2"] = betabasisk, _["rng"] = range, _["wtfd"] = wtfdjk));
                    Cmat(indexj-1, indexk-1) = Cmatjk;
                    Cmat(indexk-1, indexj-1) = Cmatjk.t();
                }
            }
            double lambdaj = betafdParj["lambda"];
            if(lambdaj > 0)
            {
               Nullable<NumericMatrix> Rmatj_ = betafdParj["penmat"];
               arma::mat Rmatj;
               if(Rmatj_.isNull()){
                List Lfdj = betafdParj["Lfd"];
                Rmatj = as<arma::mat>(eval_penalty(betabasisj, Lfdj));
               }else{
                Rmatj = as<arma::mat>(Rmatj_.get());
               }
               Cmat(indexj-1, indexj-1) = Cmat(indexj-1, indexj-1) + lambdaj * Rmatj;
            } 
        }
    }
    //ensure symmetry
    Cmat = (Cmat + Cmat.t())/2;
    arma::vec betacoef = arma::solve(Cmat, Dmat);
    //set up fdPar objects
    List betaestlist = betalist;
    mj2 = 0;
    mj1 = 0;
    for(size_t j = 0; j < p; j++)
    {
        List betafdParj = betalist[j];
        List res;
        if(betafdParj["estimate"]){
            List betafdj = betafdParj["fd"];
            List betabasisj = betafdj["basis"];
            int ncoefj = betabasisj["nbasis"];
            mj1 = mj2 + 1;
            mj2 += ncoefj;
            arma::uvec indexj = arma::conv_to<arma::uvec>::from(arma::regspace(mj1,mj2));
            arma::vec coefj = betacoef(indexj-1);
            arma::mat coefj_mat = arma::reshape(coefj, coefj.size(), 1);
            res = fd(coefj_mat, betabasisj, betafdj["fdnames"]);
            betafdParj = fdPar(res, _["lambda"]=betafdParj["lambda"]);
        }
        betaestlist[j] = betafdParj;
    }
    //set up fd objects for predicted values in yhatfdobj
    int nfine = std::max(501, 10 * nbasis + 1);
    arma::vec tfine = arma::linspace(range[0], range[1], nfine);
    arma::mat yhatmat(nfine, N, arma::fill::zeros);
    for(size_t j = 0; j < p; j++)
    {
        List xfdj = xfdlist[j];
        arma::mat xmatj = as<arma::mat>(eval_fd(_["evalarg"] = tfine, _["fdobj"] = xfdj));
        List betafdParj = betaestlist[j];
        List betafdj = betafdParj["fd"];
        arma::mat betavecj = as<arma::mat>(eval_fd(_["evalarg"] = tfine, _["fdobj"] = betafdj));
        yhatmat = yhatmat + xmatj.each_col() % arma::vectorise(betavecj); //non so se è piu efficace ma cosi va
        
    } 
    List res = smooth_basis(_["argvals"] = tfine, _["y"] = yhatmat, _["fdParobj"] = ybasisobj);
    List yhatfdobj = res["fd"];
    //betaestlist
    //yfdobj
    //
    return std::make_tuple(yfdobj, xfdlist, betalist, betaestlist, yhatfdobj, Cmat, Dmat);
}