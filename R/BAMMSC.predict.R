
BAMMSC.predict<-function(X,K,p=0.5,seed=1,...){
  
  message('### Speed up clustering process by training the BAMMSC model using portion of all cells and predict cell labels of rest cells\n')
  set.seed(seed)
  L=length(X)
  X1=X2=list()
  id1s=id2s=list()
  for(l in 1:L){
    ncell=ncol(X[[l]])
    ngene=nrow(X[[l]])
    id1=seq(ncell*p)
    id1=sample(ncell,ncell*p)
    id1=sort(id1)
    id2=setdiff(seq(ncell),id1)
    id1s[[l]]=id1
    id2s[[l]]=id2
    X1[[l]]=X[[l]][,id1]
    X2[[l]]=X[[l]][,id2]
  }
  message('### Train BAMMSC model......')
  result=BAMMSC(X1,K=K,...)
  mem1=result$mem
  
  nIter=length(result$loglik)
  alpha=Reduce('+',result$alpha[ceiling(0.9*nIter):nIter])/(nIter-ceiling(0.9*nIter)+1)
  
  message('### Predict cell cluster label......')
  mem2=list()
  for(l in 1:L){
    pos.mat=matrix(NA,ncol(X2[[l]]),K)
    for(j in 1:ncol(X2[[l]])){
      nIter=length(result$alpha)
      tmp=sapply(1:K, function(k) log(result$pi[k]) +
                   lgamma(sum(alpha[,k,l])) - lgamma(sum(X[[l]][,j])+sum(alpha[,k,l])) +
                   sum(lgamma(X2[[l]][,j] + alpha[,k,l]) - lgamma(alpha[,k,l])))
      
      for(k in 1:K){
        pos.mat[j,k] <- 1/sum(exp(tmp-tmp[k]))
      }
    }
    z=apply(pos.mat,1,which.max)
    mem2[[l]]=z
  }
  
  
  mem=list()
  for(l in 1:L){
    mem[[l]]=numeric(nrow(X[[l]]))
    mem[[l]][id1s[[l]]]=mem1[[l]]
    mem[[l]][id2s[[l]]]=mem2[[l]]
  }
  
  mem
  
}




