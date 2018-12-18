#include "RcppArmadillo.h"


// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;


arma::vec lgammaVec(arma::vec inVec){
  int G=inVec.n_elem;
  arma::vec outVec(G);
  for(int g=0; g<G; g++){
    outVec(g) = std::lgamma(inVec(g));
  }
  return outVec;
}

double ldgamma(double x,double shape,double scale){
  return (shape-1)*std::log(x)-x/scale-std::lgamma(shape)-shape*std::log(scale);
}

double sample(int K,arma::vec w){
  arma::vec x(K);
  for(int i=1;i<=K;i++) x(i-1)=i;
  arma::uvec ind=arma::sort_index(w,"descend") ;
  x=x(ind);
  w=w(ind);
  arma::vec wsum=arma::cumsum(w);
  //arma::vec a=arma::randu(1);
  int id=0;
  for(int j=0;j<w.n_elem;j++){
    //if(a(0)<=wsum(j)){
    if(Rcpp::as<double>(Rcpp::runif(1))<=wsum(j)){
      id=j;
      break;
    }
  }
  return x(id);
}



double likelihood(arma::field<arma::mat> & Xs,arma::cube & alphas,arma::field<arma::vec> & Zs,arma::field<arma::vec> & Ts,
                  int K, int G, int L, arma::mat &sigma2, arma::mat & mu,arma::vec & a, arma::vec &b){
  
  double like=0,like1=0,like2=0,like3=0,like4=0,lgamr,lgaml;
  arma::uvec ind;
  arma::vec tmp;
  
  for(int l=0;l<L;l++){
    for(int k=0;k<K;k++){
      ind =arma::find( Zs(l)==(k+1) );
      tmp=Ts(l)(ind)+arma::sum(alphas.slice(l).col(k));
      lgamr=arma::sum( std::lgamma(arma::sum(alphas.slice(l).col(k))) - lgammaVec(tmp)  );
      lgaml=0;
      for(int i=0;i<G;i++){
        tmp=Xs(l).row(i).t();
        tmp=tmp(ind)+alphas.slice(l)(i,k);
        lgaml+=arma::sum( lgammaVec( tmp ) - std::lgamma(alphas.slice(l)(i,k))  );
      }
      like1+=lgamr+lgaml;
    }
  }
  
  
  for(int k=0;k<K;k++){
    for(int i=0;i<G;i++){
      for(int l=0;l<L;l++){
        like2+= -std::log(alphas.slice(l)(i,k))-std::pow(std::log(alphas.slice(l)(i,k))-mu(i,k),2)/(2*sigma2(i,k));
      }
      like3+= -L/2*std::log(sigma2(i,k));
      like4+=ldgamma(sigma2(i,k),a(k),b(k));
    }
  }
  
  like=like1+like2+like3+like4;
  return like;
  
}



void MCMC(arma::field<arma::mat> & Xs,arma::cube & alphas,arma::field<arma::vec> & Zs,arma::field<arma::vec> & Ts,
          int K, int G, int L, arma::mat &sigma2, arma::mat & mu,arma::vec & a, arma::vec &b,
          arma::vec &Cl, arma::vec &pi_k, double sd_alpha, double sd_sigma2,
          arma::field<arma::cube> &oalphas,arma::cube &osigma2,arma::field<arma::mat> &oZs,arma::vec &olikes,arma::field<arma::cube> &odeltas,int maxIter,double likTol){
  
  double sigma2_old,sigma2_new,acc_p,lgamr_old,lgamr_new,lgaml_old,lgaml_new,alpha_new,like_new,like_old,like_sig1,like_sig2;
  arma::vec pos_p(K);
  arma::vec pos_p_temp(K);
  arma::uvec ind;
  arma::vec alpha_all_new;
  arma::vec tmp;
  int iter=0;
  double logLik = 10;
  double dif = 100.0;
  //double AIC=0;
  //double BIC=0;
  //double sumd=0;

  
  //for(int iter=0;iter<maxIter;iter++){
  while( (dif > likTol ) && (iter < maxIter) ){
    
    for(int l=0;l<L;l++){
      for(int j=0;j<Cl(l);j++){
        pos_p.fill(NA_REAL);
        pos_p_temp.fill(NA_REAL);
        for(int k=0;k<K;k++){
          pos_p_temp(k)= std::log(pi_k(k))+std::lgamma(arma::sum(alphas.slice(l).col(k)))-std::lgamma(arma::sum(Xs(l).col(j))+arma::sum(alphas.slice(l).col(k)))+
            arma::sum(lgammaVec(Xs(l).col(j)+alphas.slice(l).col(k))-lgammaVec(alphas.slice(l).col(k)));
        }
        for(int k=0;k<K;k++){
          pos_p(k)=1/arma::sum(arma::exp(pos_p_temp-pos_p_temp[k]));
          odeltas(l)(j,k,iter)=pos_p(k);
        }
        Zs(l)(j)=sample(K,pos_p);
        oZs(l)(j,iter)=Zs(l)(j);
        
      }
    }
    
    
    for(int i=0;i<G;i++){
      for(int k=0;k<K;k++){
        sigma2_old=sigma2(i,k);
        sigma2_new=Rcpp::as<double>(Rcpp::rnorm(1,sigma2_old,sd_sigma2));
        if(sigma2_new>100) sigma2_new=200-sigma2_new;
        if(sigma2_new<0) sigma2_new= -sigma2_new;
        like_sig1=0;like_sig2=0;
        for(int l=0;l<L;l++){
          like_sig1+= -std::pow(std::log(alphas.slice(l)(i,k))-mu(i,k),2)/(2*sigma2_new);
          like_sig2+= -std::pow(std::log(alphas.slice(l)(i,k))-mu(i,k),2)/(2*sigma2_old);
        }
        like_sig1=like_sig1-L/2*std::log(sigma2_new)+ldgamma(sigma2_new,a(k),b(k));
        like_sig2=like_sig2-L/2*std::log(sigma2_old)+ldgamma(sigma2_old,a(k),b(k));
        acc_p=like_sig1-like_sig2;
        if( acc_p>=0 ) acc_p=0;
        if(std::log(Rcpp::as<double>(Rcpp::runif(1)))<acc_p){
          sigma2(i,k)=sigma2_new;
        }
      }
    }
    osigma2.slice(iter)=sigma2;
    
    for(int l=0;l<L;l++){
      for(int k=0;k<K;k++){
        ind =arma::find( Zs(l)==(k+1) );
        for(int i=0;i<G;i++){
          alpha_all_new=alphas.slice(l).col(k);
          lgamr_old=arma::sum( std::lgamma(arma::sum(alphas.slice(l).col(k)))-lgammaVec(Ts(l)(ind)+arma::sum(alphas.slice(l).col(k))) ); //variant
          tmp=Xs(l).row(i).t();
          tmp=tmp(ind)+alphas.slice(l)(i,k);
          lgaml_old=arma::sum( lgammaVec(tmp)-std::lgamma(alphas.slice(l)(i,k)) ); //invariant
          alpha_new=Rcpp::as<double>(Rcpp::rnorm(1,alphas.slice(l)(i,k),sd_alpha));
          if(alpha_new>400) alpha_new=800-alpha_new;
          if(alpha_new<0) alpha_new= -alpha_new;
          alpha_all_new(i)=alpha_new;
          lgamr_new=arma::sum( std::lgamma(arma::sum(alpha_all_new))-lgammaVec(Ts(l)(ind)+arma::sum(alpha_all_new)) ); //variant
          
          tmp=Xs(l).row(i).t();
          tmp=tmp(ind)+alpha_new;
          lgaml_new=arma::sum( lgammaVec(tmp)-std::lgamma(alpha_new) ); //invariant
          like_new=lgaml_new+lgamr_new-std::log(alpha_new)-std::pow(std::log(alpha_new)-mu(i,k),2)/(2*sigma2(i,k)); //variant
          like_old=lgaml_old+lgamr_old-std::log(alphas.slice(l)(i,k))-std::pow(std::log(alphas.slice(l)(i,k))-mu(i,k),2)/(2*sigma2(i,k)); //variant
          acc_p=like_new-like_old;
          if(acc_p>=0) acc_p=0;
          if(std::log(Rcpp::as<double>(Rcpp::runif(1)))<acc_p){
            alphas.slice(l)(i,k)=alpha_new;
          }
        }
      }
    }
    oalphas(iter)=alphas;
    olikes(iter)=likelihood(Xs,alphas,Zs,Ts,K,G, L,sigma2,mu,a,b);
    
    dif = std::abs((olikes(iter) - logLik)/logLik);
    logLik=olikes(iter);
    iter++;
    
    if(iter%100==1) std::cout<<"iteration:"<<iter<<std::endl;
  }

  
  //AIC = (-2)*olikes(iter) + 2*((K*G)+(K*J));
  //BIC = (-2)*olikes(iter) + std::log(G*J)*((K*G)+(K*J));
  
}



RcppExport SEXP MCMC_multinomial(SEXP Xs_, SEXP alphas_,SEXP Zs_,SEXP Ts_,
                                 SEXP K_, SEXP G_, SEXP L_,
                                 SEXP sdalpha_, SEXP sdsigma2_, 
                                 SEXP sigma2_, SEXP mu_,  SEXP a_,  SEXP b_,  SEXP Cl_,  SEXP pik_, SEXP maxIter_, SEXP likTol_) {
  
  
  Rcpp::List Xs(Xs_);
  Rcpp::List alphas(alphas_);
  Rcpp::List Zs(Zs_);
  Rcpp::List Ts(Ts_);

  int K = Rcpp::as<int>(K_);
  int G = Rcpp::as<int>(G_);
  int L = Rcpp::as<int>(L_);
  double sdalpha=Rcpp::as<double>(sdalpha_);
  double sdsigma2=Rcpp::as<double>(sdsigma2_);
  
  arma::mat sigma2 = Rcpp::as<arma::mat>(sigma2_);
  arma::mat mu = Rcpp::as<arma::mat>(mu_);
  arma::vec a = Rcpp::as<arma::vec>(a_);
  arma::vec b = Rcpp::as<arma::vec>(b_);
  arma::vec Cl = Rcpp::as<arma::vec>(Cl_);
  arma::vec pik = Rcpp::as<arma::vec>(pik_); //dimension?
  
	int maxIter = Rcpp::as<int>(maxIter_);
  double likTol=Rcpp::as<double>(likTol_);
  
  arma::field<arma::mat> arma_Xs(L);
  arma::field<arma::vec> arma_Zs(L);
  arma::field<arma::vec> arma_Ts(L);
  arma::cube arma_alphas(G,K,L);
  arma::vec olikes(maxIter);

  for (int i=0; i<L; i++) {
    arma_Xs(i) = Rcpp::as<arma::mat>(Xs[i]);
    arma_alphas.slice(i)=Rcpp::as<arma::mat>(alphas[i]);
    arma_Zs(i) = Rcpp::as<arma::vec>(Zs[i]);
    arma_Ts(i) = Rcpp::as<arma::vec>(Ts[i]);
  }

  //output variables
  arma::field<arma::cube> oalphas(maxIter);
  for (int i=0; i<maxIter; i++) {
    oalphas(i)=arma_alphas;
    oalphas(i).fill(NA_REAL);
  }
  
  arma::field<arma::cube> odeltas(L);
  for (int i=0; i<L; i++) {
    arma::cube arma_delta(Cl(i),K,maxIter);
    odeltas(i)=arma_delta;
    odeltas(i).fill(NA_REAL);
  }
  
  
  arma::cube osigma2(G,K,maxIter);
  osigma2.fill(NA_REAL);
  arma::field<arma::mat> oZs(L);

  for (int i=0; i<L; i++) {
    oZs(i)=arma::mat(arma_Zs(i).n_elem,maxIter);
    oZs(i).fill(NA_REAL);
  }

  MCMC(arma_Xs,arma_alphas,arma_Zs,arma_Ts,K,G,L,sigma2,mu,a,b,Cl,pik,sdalpha,sdsigma2,oalphas,osigma2,oZs,olikes,odeltas,maxIter,likTol);
   
  Rcpp::List result;
  result["mem"]=oZs;
  result["loglik"]=olikes;
  result["sigma2"]=osigma2;
  result["alpha"]=oalphas;
  result["pi"]=pik;
  result["delta"]=odeltas;
  //result["AIC"] = AIC;
  //result["BIC"] = BIC;
  
  return result;

}







  
