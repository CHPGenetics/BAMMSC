The following R scripts are used to reproduce the results for Human PBMC dataset used in the manuscript. 
The processed Human PBMC data could be downloaded from the link [HumanPBMC](https://drive.google.com/open?id=1t0li4AalSRNpdLnGljhj5wJFB0JPlvm7)

```
# load packages
library(mclust)
library(BAMMSC)
library(cellrangerRkit)


# load datasets
load("Human_PBMC.rda")
load("HumanPBMC_ApproxTruth.Rdata")
load("HumanPBMC_PooledTSNE.Rdata")

# run BAMMSC
set.seed(12352)
result<-BAMMSC(count, K=7,nBurn = 500, maxIter.BAMMSC = 1500)
adjustedRandIndex(unlist(result$mem)[ApproxTruth[[1]]],ApproxTruth[[2]])

# run t-SNE, colored by clustering label
clust<-unlist(result$mem)
tsne_clust <- data.frame(Barcode=tsne_result[,1],TSNE.1=tsne_result[,2],TSNE.2=tsne_result[,3],Clust=clust)

pdf(file="HumanPBMC_BAMMSC.pdf")
visualize_clusters(as.character(tsne_clust$Clust),tsne_clust[c("TSNE.1","TSNE.2")],title="BAMM-SC clustering",marker_size=0.5)
dev.off()

# run t-SNE, colored by Sample ID #
batch<-c(rep("Sample 1",dim(count[[1]])[2]),rep("Sample 2",dim(count[[2]])[2]),
         rep("Sample 3",dim(count[[3]])[2]),rep("Sample 4",dim(count[[4]])[2]),
         rep("Sample 5",dim(count[[5]])[2]))
tsne_clust <- data.frame(Barcode=tsne_result[,1],TSNE.1=tsne_result[,2],
                         TSNE.2=tsne_result[,3],Batch=batch)
pdf(file="HumanPBMC_Sample.pdf")
visualize_clusters(as.character(tsne_clust$Batch),tsne_clust[c("TSNE.1","TSNE.2")],
                   title="Sample Labels",marker_size=0.5)
dev.off()

```

<img src="figures/HumanPBMC_BAMMSC.png" style="display: block; margin: auto;" />




