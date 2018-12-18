

BAMMSC <- function(X, K=4,option='BAMMSC',method_cluster_initial="kmeans", method_alpha_initial="Ronning",maxIter.DIMMSC=400,tol.DIMMSC=1e-4, likTol.DIMMSC=1e-2, customized_initial=NULL,
                   nBurn=400,maxIter.BAMMSC=500,likTol.BAMMSC=1e-10,trace=F){

  
  message('Calculate the initial values......')
  
  low.bound=0.001
  sdalpha=1
  sdsigma2=0.5
  
  if(class(X)=='data.frame' || class(X)=='matrix'){
    message('Only one individual, call DIMMSC......')
    res=DIMMSC(X,K,method_cluster_initial=method_cluster_initial,method_alpha_initial=method_alpha_initial,maxIter=maxIter.DIMMSC, tol=tol.DIMMSC, lik.tol=likTol.DIMMSC,customized_initial=customized_initial,trace=trace)
    return(list(mem=res$mem,loglik=res$loglik,AIC=res$AIC,BIC=res$BIC,alpha=res$alpha,pi=res$pi,delta=res$delta))
  }else if(class(X)=='list'){
    message('Multiple individuals, call BAMMSC......')
    L=length(X)
    Ts=list()
    for (l in 1:L){
      Ts[[l]]<-colSums(X[[l]])
    }
    G=nrow(X[[1]])
    Cl=sapply(X,ncol)
    
    XX=do.call(cbind,X)
    res.all=DIMMSC(XX,K,method_cluster_initial=method_cluster_initial,method_alpha_initial=method_alpha_initial,maxIter=maxIter.DIMMSC, tol=tol.DIMMSC, lik.tol=likTol.DIMMSC,customized_initial=customized_initial,trace=trace)
    if(option=='DIMMSC'){
      message('Multiple individuals merged, call DIMMSC......')
      return(list(mem=res.all$mem,loglik=res.all$loglik,AIC=res.all$AIC,BIC=res.all$BIC,alpha=res.all$alpha,pi=res.all$pi,delta=res.all$delta))
    }
    
    res.all$alpha[res.all$alpha<low.bound]=low.bound
    
    alphas=list()
    for (l in 1:L){
      alphas[[l]]<-t(res.all$alpha)
    }
    
    Zs=list()
    cc=c(0,cumsum(Cl))
    for(l in 1:L){
      Zs[[l]]<-res.all$mem[(1+cc[l]):cc[l+1]]
    }
    
    alpha.L<-list()
    for( l in 1:L){
      res=DIMMSC(X[[l]],K,method_cluster_initial=method_cluster_initial, method_alpha_initial=method_alpha_initial,maxIter=maxIter.DIMMSC,tol=tol.DIMMSC, lik.tol=likTol.DIMMSC,customized_initial=customized_initial,trace=trace)
      alpha.L[[l]]<-t(res$alpha)
      alpha.L[[l]][alpha.L[[l]]<low.bound]=low.bound
    }
    
    
    match.id=match_cell(alpha=alpha.L,L=L,K=K)
    alpha.L.order=list()
    for(l in 1:L){
      alpha.L.order[[l]]=alpha.L[[l]][,match.id[,l]]
    }
    
    
    sigma2=matrix(0,G,K)
    mu=matrix(0,G,K)
    
    for(g in 1:G){
      for(k in 1:K){
        sig=numeric(L)
        for(l in 1:L){
          sig[l]=alpha.L.order[[l]][g,k]
        }
        if(var(sig)==0){
          sigma2[g,k]=low.bound
        }else{
          sigma2[g,k]=var(log(sig))
        }
        mu[g,k]=mean(log(sig))
      }
      if(var(sigma2[,k])==0){
        sigma2[,k]<-rgamma(G,shape=0.1,scale=0.1)
      }
      sigma2[sigma2[,k]<low.bound,k]<-low.bound
    }
    
    
    a=b=matrix(0,1,K)
    for(k in 1:K){
      b[k]=var(sigma2[,k])/mean(sigma2[,k])
      a[k]=mean(sigma2[,k])/b[,k]
    }
    
    message('MCMC performs.....')
    res=MCMC_multinomial(X, alphas,Zs, Ts,K,G,L,sdalpha,sdsigma2,sigma2, mu, a, b, Cl, res.all$pi, maxIter.BAMMSC, likTol.BAMMSC)
    
    mems=list()
    for(l in 1:L){
      mems[[l]]=numeric(Cl[l])
      for(j in 1:Cl[l]){
        mems[[l]][j]=Mode(res$mem[[l]][j,nBurn:maxIter.BAMMSC])
      }
    }
    
    loglik.burnin=mean(res$loglik[nBurn:maxIter.BAMMSC])
    AIC=(-2)*loglik.burnin+2*(sum(Cl)*K+K*G*L+K*G)
    BIC=(-2)*loglik.burnin+log(sum(Cl)*G)*(sum(Cl)*K+K*G*L+K*G)
    
    return (list(mem=mems,loglik=res$loglik,AIC=AIC,BIC=BIC,alpha=res$alpha,pi=res$pi,sigma2=res$sigma2,delta=res$delta))
    
  }
  
}









