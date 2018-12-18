#include <RcppArmadillo.h>


// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


arma::vec lgammaVec2(arma::vec inVec){
  int G=inVec.n_elem;
  arma::vec outVec(G);
  for(int g=0; g<G; g++){
    outVec(g) = std::lgamma(inVec(g));
     if(!isfinite(outVec(g))){ // if inVec(g) is closer than -1e-154 to 0, e.g. -1e-155, -inf will be returned.
            //So we would truncate it to a certain number: here it's lgamma(-1e-154)
            outVec(g) = 354.5891; // lgamma(-1e-154)
     }
  }
  return outVec;
}


arma::vec rdirichlet (int n){
    arma::vec v = arma::randg<arma::vec>(n);
    double sum = arma::sum(v);
    arma::vec out = v/sum;
    return(out);
}

int indexMax( arma::vec& inVec){
    int n = inVec.size();
    int idx = 0;
    double tmpValue = inVec(0);
    if(n>1){
        for(int i=1; i<n; i++){
            if(inVec(i) > tmpValue){
                idx = i;
                tmpValue = inVec(i);
            }
        }
    }
    return(idx);
}



// [[Rcpp::export]]
RcppExport SEXP EM_multinomial(SEXP X_, SEXP K_, SEXP alpha_, SEXP maxIter_, SEXP tol_, SEXP likTol_){


  arma::mat X = Rcpp::as<arma::mat>(X_);
  int K = Rcpp::as<int>(K_);
  arma::mat alphaInput = Rcpp::as<arma::mat>(alpha_);
  int maxIter = Rcpp::as<int>(maxIter_);
  double tol=Rcpp::as<double>(tol_);
  double likTol=Rcpp::as<double>(likTol_);
//Rcpp::List EM_multinomial(const arma::mat& X, const int K, const arma::mat& alphaInput, const int maxIter, const double tol, const double likTol){   //X[G,J], alpha[K,G]
    arma::mat alpha=alphaInput;
    arma::vec repTmp(K);
    repTmp.fill(1);
    arma::vec pi(K);
    pi = rdirichlet(K);
    arma::vec piNew(K);

    int J = X.n_cols;
    int G = X.n_rows;
    double differ = 1.0;
    int iter = 0;
    double logLik = 1.0;
    double dif = 100.0;
    arma::mat delta = arma::zeros<arma::mat>(J,K);

    while( (differ > tol || dif > likTol) && (iter < maxIter) ){
        // E-step compute omega
        arma::mat num1 = arma::zeros<arma::mat>(J,K);
        arma::mat num2 = arma::zeros<arma::mat>(J,K);

        arma::mat lgammaAlpha = arma::zeros<arma::mat>(G,K); // lgamma(alpha[,].t())
        arma::vec lgammaSumAlpha(K); //lgamma(sum(alpha(k,)))
        arma::vec dataAlphaTmp(G); //X[,j]+alpha[,k]
        arma::vec logpi(K); //log(pi(k))
        arma::mat deltaTmp = arma::zeros<arma::mat>(J,K);
        arma::mat alphaT = alpha.t();
        arma::vec alphaSum(K);
        arma::vec dataSum(J);
        delta.fill(0);

        for(int k=0; k<K; k++){
            // arma::vec tmpV = alpha.row(k).t();//always remember to transpose rowvec when you are surely using vec.
            arma::vec tmpV = alphaT.col(k);
            lgammaAlpha.col(k) = lgammaVec2(tmpV);
            lgammaSumAlpha(k) = lgamma(arma::sum(alpha.row(k)));
            logpi(k) = log(pi(k));
            alphaSum(k) = arma::sum(alpha.row(k));
        }
        for(int j=0; j<J; j++){
            dataSum(j) = arma::sum(X.col(j));
            for(int k=0; k<K; k++){
                dataAlphaTmp = X.col(j) + alphaT.col(k);
                num1(j,k) = arma::sum( lgammaVec2(dataAlphaTmp) - lgammaAlpha.col(k) );
                num2(j,k) = lgammaSumAlpha(k) - lgamma( alphaSum(k) + dataSum(j) );
                deltaTmp(j,k) = num1(j,k) + num2(j,k) + logpi(k);
            }
        }

            for(int j=0; j<J; j++){
                for(int k=0; k<K; k++){
                    double sumTmp=0.0;
                    for(int m=0; m<K; m++){
                        if(m!=k){
                            sumTmp += exp(arma::as_scalar(deltaTmp(j,m) - deltaTmp(j,k)));
                        }
                    }
                    delta(j,k) = 1.0/(1.0+sumTmp);
                }
            }

            //M-step: update pi and alpha

            for(int k=0; k<K; k++){
                piNew(k) = arma::sum(delta.col(k))/double(J);
            }

            arma::mat alphaNew = arma::zeros<arma::mat>(K,G);

            for(int k=0; k<K; k++){
                double den = 0.0;
                for(int j=0; j<J; j++){;
                    den += arma::as_scalar(delta(j,k))*dataSum(j)/( dataSum(j)-1+alphaSum(k) );
                }
                for(int g=0; g<G; g++){
                    double tmpNum=0.0;
                    arma::vec dataG(J);
                    dataG = ((X.row(g)+0.000001) / (X.row(g)+0.000001-1+alpha(k,g))).t();
                    // tmpNum += arma::as_scalar(delta.col(k).t()*dataG); //after matrix operation it stays a matrix
                    // below is another option of above, fix, are the results same?
                    for(int j=0; j<J; j++){
                        tmpNum += delta(j,k)*dataG(j);
                    }
                    alphaNew(k,g) = alpha(k,g)*tmpNum/den;
                }
            }

            arma::uvec mem(J);
            for(int j=0; j<J; j++){
                arma::vec tmpVec = delta.row(j).t();
                mem(j) = indexMax(tmpVec); //.index_max(); //elements are unsigned int
            }
            // arma::uvec sort = arma::sort_index(pi);
            arma::vec sort(K); //actually is rank function in R
            arma::uvec sortidx = arma::sort_index(pi);
            int tmpk=0;
            while(tmpk<K){
                for(int k=0; k<K; k++){
                    if(sortidx(k)==tmpk){
                        sort(tmpk)=k;
                        tmpk++;
                        break;
                    }
                }
            }
            arma::uvec res = mem;
            for(int k=0; k<K; k++){
                for(int j=0; j<J; j++){
                    if(mem(j)==k){
                        res(j)=sort(k);
                    }
                }
            }
            mem=res;
            arma::mat num=num1+num2;
            arma::vec lik(J);
            for(int j=0; j<J; j++){
                lik(j) = num(j,mem(j));
            }
            double newLogLik = arma::sum(lik);

            // calculate diff to check convergence
            dif = std::abs((newLogLik - logLik)/logLik*100.0);

            double sumd=0.0;
            for(int k=0; k<K; k++){
                sumd += pow((arma::as_scalar(piNew(k) - pi(k))), 2.0);
            }

            differ = sqrt(std::abs(sumd));
            pi=piNew;
            alpha=alphaNew;
            logLik = newLogLik;

            for(int k=0; k<K; k++){
                for(int g=0; g<G; g++){
                    if(alpha(k,g)==0){
                        alpha(k,g) = 0.000001;
                    }
                }
            }

            iter++;
    }

        arma::uvec mem(J);
        for(int j=0; j<J; j++){
            for(int k=0; k<K; k++){
                arma::vec tmpVec = delta.row(j).t();
                mem(j) = indexMax(tmpVec); //.index_max() will use the last index when max is tie.
            }
        }
        double AIC, BIC;
        AIC = (-2)*logLik + 2*((K*G)+(K*J));
        BIC = (-2)*logLik + log(G*J)*((K*G)+(K*J));

        Rcpp::List result;

        result["pi"] = pi;
        result["delta"] = delta;
        result["alpha"] = alpha;
        result["mem"] = mem + 1;
        result["loglik"] = logLik;
        result["AIC"] = AIC;
        result["BIC"] = BIC;

        return result;
}

