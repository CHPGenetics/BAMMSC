# BAMMSC
A Bayesian mixture model for clustering droplet-based single cell transcriptomic data from population studies

## Description
BAMMSC is an R package for clustering droplet-based single cell transcriptomic data from multiple individuals simultaneously. It adopts a Bayesian hierarchical Dirichlet multinomial mixture model, which explicitly characterizes three levels of variabilities (i.e., genes, cell types and individuals). BAMMSC is able to take account for data heterogeneity and batch effect, such as unbalanced sequencing depths, variable read length and hidden technical bias, among multiple individuals. BAMMSC also integrates DIMMSC for single individual analysis.

# Author
Zhe Sun, Li Chen, Qianhui Huang, Anthony Richard Cillo, Tracy Tabib, Ying Ding, Jay Kolls, Robert Lafyatis, Dario Vignali, Kong Chen, Ming Hu,* and Wei Chen* [link](https://www.biorxiv.org/content/biorxiv/early/2018/08/16/392662.full.pdf)

# Maintainer
Zhe Sun <zhs31@pitt.edu>, Li Chen <li.chen@auburn.edu>


# Install BAMMSC
```r
install.packages("devtools")
library(devtools)
install_github("lichen-lab/BAMMSC")
```


# Descriptions for BAMMSC

## Usage
BAMMSC(X, K=4, option="BAMMSC",method_cluster_initial="kmeans", method_alpha_initial="Ronning", maxIter.DIMMSC=400, tol.DIMMSC=1e-4, likTol.DIMMSC=1e-2, nBurn=400, maxIter.BAMMSC=500, likTol.BAMMSC=1e-10)

## Arguments
*  X: a list with each element a UMI count matrix for each individual, with row as number of genes and column as number of cells. The numbers of genes are the same across different individuals but the number of cells might differ. If X is a matrix, DIMMSC will be used
*  K: number of clusters, default is 4
*  Option: Default is  "BAMMSC"; if "DIMMSC", all cells will be combined across individuals and DIMMSC will be used
*  method_cluster_initial: method for initializing clusters in the procedure of DIMMSC, "kmeans" (default) or "random"
*  method_alpha_initial: method for initializing the alpha matrix in the procedure of DIMMSC, "Ronning" (default, Ronning's method, 1989) or "Weir" (Weir and Hill's method, 2002)
*  maxiter.DIMMSC: maximum number of iterations for DIMMSC, default is 400
*  tol.DIMMSC: a convergence tolerance for the difference of vector pie between iterations in the procedure of DIMMSC, default is 1e-2
*  likTol.DIMMSC: a convergence tolerance for the difference of log-likelihoods between iterations in the procedure of DIMMSC, default is 1e-4
*  nBurn: number of iterations discarded at the initial portion of the Markov chain, default is 400
*  maxIter.BAMMSC: maximum number of iterations in MCMC, default is 500
*  likTol.BAMMSC: a convergence tolerance for the difference of log-likelihoods between iterations in MCMC, default is 1e-10

## Output values
* mem: a list of vectors. Each vector [#cell] is the clustering membership for all cells in each individual. Length of list is #individuals
* loglik: a vector of log likelihood all iteration in MCMC
* AIC: Akaike information criterion (AIC)
* BIC: Bayesian information criterion (BIC)
* alpha: a list of 3-D matrix. Each 3-D matrix [#gene, #cluster, #individuals] is alpha estimates for each iteration. Length of list is #iterations
* pi: a vector of pi estimates
* sigma2: a 3-D matrix of sigma^2 estimates [#gene, #cluster, #iteration]
* delta: a list of 3-D matrix. Each matrix [#cell, #cluster, #iteration] is the poster probability estimates for each individual across all iterations. Length of list is #individuals

# Descriptions for BAMMSC.predict

## Usage
BAMMSC.predict(X,K=4,p=0.5,...)

## Arguments
*  X: a list with each element a UMI count matrix for each individual, with row as number of genes and column as number of cells. The numbers of genes are the same across different individuals but the number of cells might differ. If X is a matrix, DIMMSC will be used
*  K: number of clusters, default is 4
*  p: proportion of cells used to train a BAMMSC model and the cluster labels of the rest cells are predicted by the trained label
 *  ...: other arguments the same as BAMMSC.
 
 ## Output values
 * mem: a list of vectors. Each vector [#cell] is the clustering membership for all cells in each individual. Length of list is #individuals


# Examples

Use BAMMSC and DIMMSC
```r
# Load the example data data_BAMMSC
# data_BAMMSC contains 10 individuals with 100 genes and 400 cells each
library(BAMMSC)
data("data_BAMMSC")
# run BAMMSC across multiple individuals
result=BAMMSC(data_BAMMSC,K=4)

# Load the example data data_DIMMSC
# data_DIMMSC contains 1000 genes and 6000 cells.
data("data_DIMMSC")
# For single individual, BAMMSC could call DIMMSC in two ways.
result=DIMMSC(data_DIMMSC,K=3)
result=BAMMSC (data_DIMMSC,K=3)
```

Use BAMMSC.predict and compare to BAMMSC
```r
library(BAMMSC)
library(mclust)
data("data_BAMMSC")

#create true cell labels for all cells in 10 individuals
K=4
z=rep(seq(K),each=100)
Z=list(z,z,z,z,z,z,z,z,z,z)  

#obtain all cell labels using BAMMSC.predict
Zpred=BAMMSC.predict(X=data_BAMMSC,K=4,p=0.5,seed=1)
adjustedRandIndex(unlist(Zpred),unlist(Z))

#obtain all cell labels using BAMMSC
result=BAMMSC(X=data_BAMMSC,K=4)
adjustedRandIndex(unlist(result$mem),unlist(Z))

```

# Use BAMMSC for 10X Genomics 
Tutorial to use BAMMSC for 10X Genomics could be found at https://github.com/lichen-lab/BAMMSC/blob/master/10XGenomics.md


# Use BAMMSC to run three real datasets

The Human PBMC data analysis could be found at https://github.com/lichen-lab/BAMMSC/blob/master/HumanPBMC.md

The Human skin data analysis could be found at https://github.com/lichen-lab/BAMMSC/blob/master/Humanskin.md

The Mouse lung data analysis could be found at https://github.com/lichen-lab/BAMMSC/blob/master/Mouselung.md

# Fastq files of three real datasets

All fastq files could be downloaded [here](https://pitt.app.box.com/s/aqwl3aedfqxp41oecvudybv9h9u0t4ai)





