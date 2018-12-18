
DIMMSC <- function(X, K=2, method_cluster_initial="kmeans", method_alpha_initial="Ronning", maxIter=200, tol=1e-4, lik.tol=1e-2, customized_initial=NULL,trace=TRUE)
{
  # Initialize clusters
  if(trace) cat('Initializing clusters ')
  if(method_cluster_initial == "kmeans"){
    if(trace) cat('with kmeans clusters.\n')
    kX <- as.matrix(log2(X+1)) # Normalize X for k-means
    k <- kmeans(t(kX),K)
    clusters_initial <- k$cluster
  } else if(method_cluster_initial == "random"){
    if(trace) cat('with random clusters.\n')
    clusters_initial <- sample(seq(K), ncol(X), replace = T)
  } else if(method_cluster_initial == "customized"){
    if(trace) cat('with customized clusters.\n')
    if(!(length(customized_initial) == ncol(X) & max(as.numeric(names(table(customized_initial)))) == K & min(as.numeric(names(table(customized_initial)))) == 1)){
      stop('Customized initial clusters should be a vector of positive integers {1,2,...,K}, the length should be equal to the number of cells.')
    } else if(!(sum(as.numeric(names(table(customized_initial))) == seq(K)) == K)){
      stop('Customized initial clusters should be a vector of positive integers {1,2,...,K}, the length should be equal to the number of cells.')
    } else{
      clusters_initial <- customized_initial
    }
  } else{
    stop('Method for initializing clusters should be "kmeans", "random" or "customized".')
  }
  
  # Initialize alpha matrix
  if(!(method_alpha_initial == "Ronning" | method_alpha_initial == "Weir")){
    stop('Method for initializing alpha matrix should be "Ronning" or "Weir".')
  }
  if(trace) cat('Initializing alpha matrix with',method_alpha_initial,'method.\n')
  alpha <- EM_initial_alpha(X=X, clusters_initial=clusters_initial, method_alpha_initial=method_alpha_initial)
  
  # DIMMSC clustering
  if(trace) cat('Performing DIMMSC clustering...\n')
  result <- EM_multinomial(X=X, K=K, alpha=alpha, maxIter=maxIter, tol=tol, likTol=lik.tol)
  
  if(trace) cat('Analysis is finished.\n')
  result
}










