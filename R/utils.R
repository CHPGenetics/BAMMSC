EM_initial_alpha <- function(X, clusters_initial, method_alpha_initial)
{
  num_cluster <- length(unique(clusters_initial))
  cluster <- matrix(,num_cluster,1)
  for (i in 1:num_cluster){
    cluster[i,] <-  length(which(clusters_initial==i))/ncol(X)
  }
  cluster <- as.numeric(cluster)
  sort <- rank(cluster,ties.method="random")
  label.new <- clusters_initial
  for(j in 1:num_cluster){
    label.new[which(clusters_initial==j)] <- sort[j]
  }
  
  location <- list()
  for (m in 1:num_cluster){
    location[[m]] <- which(label.new==m)
  }
  newX <- list()
  for (m in 1:num_cluster){
    newX[[m]] <- X[,location[[m]]]
  }
  
  p <- matrix(,nrow(X),num_cluster)
  for (m in 1:num_cluster){
    sum <- sum(as.vector(newX[[m]]))
    for ( i in 1:nrow(X)){
      p[i,m] <- sum(newX[[m]][i,])/sum
    }
  }
  
  new.X <- list()
  for (m in 1:num_cluster){
    new.X[[m]] <- newX[[m]][rowSums(newX[[m]]) != 0, colSums(newX[[m]]) != 0]
  }
  C <- matrix(,num_cluster,1)
  for (m in 1:num_cluster){
    C[m,1] <- ncol(new.X[[m]])
  }
  G <- matrix(,num_cluster,1)
  for (m in 1:num_cluster){
    G[m,1] <- nrow(new.X[[m]])
  }
  
  ppp <- list()
  for (m in 1:num_cluster){
    ppp[[m]] <- matrix(,G[m,1],C[m,1])
    for (j in 1:C[m,1]){
      tmp = sum(new.X[[m]][,j])
      ppp[[m]][,j]=new.X[[m]][,j]/tmp
    }
  }
  
  pp <- list()
  for (m in 1:num_cluster){
    pp[[m]] <- matrix(,G[m,1],1)
    for ( i in 1:G[m,1]){
      pp[[m]][i,1] <- mean(ppp[[m]][i,])
    }
  }
  
  v <- list()
  for (m in 1:num_cluster){
    v[[m]] <- matrix(,G[m,1],1)
    for ( i in 1:G[m,1]){
      v[[m]][i,1] <- var(ppp[[m]][i,])
    }
    v[[m]][which(v[[m]]==0),1] <- mean(v[[m]],na.rm=T)
  }
  
  s <- matrix(,num_cluster,1)
  for(m in 1:num_cluster){
    sum <- 0
    for (i in 1:(G[m,1]-1)){
      tmp <- log( ( pp[[m]][i,1]*(1-pp[[m]][i,1])/ v[[m]][i,1] ) -1 )
      sum <- sum+tmp
    }
    s[m,1] <- exp( (1/(G[m,1]-1))*sum )
  }
  
  new.alpha <- list()
  for(m in 1:num_cluster){
    new.alpha[[m]] <- s[m,1]*p[,m]
    new.alpha[[m]][which(new.alpha[[m]]==0)] <- 0.000001
  }
  
  new_alpha_R <- matrix(,num_cluster,nrow(X))
  for (m in 1:num_cluster){
    new_alpha_R[m,] <- new.alpha[[m]]
  }
  
  new.alpha <- list()
  for (m in 1:num_cluster){
    mom <- weirMoM(t(new.X[[m]]), se=FALSE)
    if (mom <= 0) {mom <- 0.005}
    initscalar <- (1 - mom)/mom
    new.alpha[[m]] <- p[,m]*initscalar
    new.alpha[[m]][which(new.alpha[[m]]==0)] <- 0.000001
  }
  
  new_alpha_W <- matrix(,num_cluster,nrow(X))
  for (m in 1:num_cluster){
    new_alpha_W[m,] <- new.alpha[[m]]
  }
  
  if(method_alpha_initial == "Ronning"){
    return(new_alpha_R)
  } else{
    return(new_alpha_W)
  }
}


match_cell<-function(alpha,L,K){
  
  cor<-list()
  for(l in 1:L) cor[[l]]<-matrix(,K,K)
  for(l in 1:L){
    for(i in 1:K){
      for(j in 1:K){
        cor[[l]][i,j]<-cor(alpha[[1]][,i],alpha[[l]][,j])
      }
    }
  }
  
  l1<-list()
  for(l in 1:L) l1[[l]]<-matrix(,K,K)
  for(l in 1:L){
    for(i in 1:K){
      for(j in 1:K){
        l1[[l]][i,j]<-norm(as.matrix(alpha[[1]][,i])-as.matrix(alpha[[l]][,j]),type="1")
      }
    }
  }
  
  max_l1<-matrix(,K,L)
  for(i in 1:K){
    for(j in 1:L){
      max_l1[i,j]<-which(l1[[j]][i,]==min(l1[[j]][i,]))
    }
  }
  
  check<-c()
  for(i in 1:L){
    left<-max_l1[,1][is.na(match(max_l1[,1],max_l1[,i]))]
    check<-c(check,left)
  }
  
  while(sd(as.numeric(max_l1))!=sd(rep(c(1:K),L))|(length(check)!=0 )) {
    check<-c()
    for(i in 1:L){
      left<-max_l1[,1][is.na(match(max_l1[,1],max_l1[,i]))]
      dup<- unique(max_l1[,i][duplicated(max_l1[,i])])
      select<-which(duplicated(max_l1[,i])|duplicated(max_l1[,i][K:1])[K:1])
      
      if (length(left)>1){
        for(j in 1:length(left)){
          sort<-order(cor[[i]][select,left[j]])
          sel<-select[sort[-1]]
          max_l1[sel,i]<-left[1]
          select<-sel
        }
      }
      
      if(length(left)==1){
        sort<-order(cor[[i]][select,left])
        sel<-select[sort[-1]]
        max_l1[sel,i]<-left
        select<-sel
      }
    }
  }
  max_l1
}


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}



MCMC_multinomial <- function(Xs, alphas, Zs, Ts, K, G,L,sdalpha,sdsigma2,sigma2, mu, a, b, Cl, pik, maxIter, likTol) {
  .Call('MCMC_multinomial', PACKAGE = 'BAMMSC', Xs, alphas, Zs, Ts, K, G,L,sdalpha,sdsigma2,sigma2, mu, a, b, Cl, pik, maxIter, likTol)
}


EM_multinomial <- function(X, K, alpha, maxIter, tol, likTol) {
  .Call('EM_multinomial', PACKAGE = 'BAMMSC', X, K, alpha, maxIter, tol, likTol)
}












