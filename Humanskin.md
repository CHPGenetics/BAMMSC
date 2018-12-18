The following R scripts are used to reproduce the results for Human skin dataset used in the manuscript. 
The processed Human skin data could be downloaded from the link [Humanskin](https://drive.google.com/open?id=1Dk4CLqSdkoR5E00zzV0wETlFZO1lCKLv)

```
# load packages
library(mclust)
library(BAMMSC)
library(cellrangerRkit)

# load datasets
load("Human_skin.rda")
load("HumanSkin_ApproxTruth.Rdata")
load("HumanSkin_PooledTSNE.Rdata")

# run BAMMSC
set.seed(14)
result<-BAMMSC(count ,K=8,nBurn = 800, maxIter.BAMMSC = 1500)
adjustedRandIndex(unlist(result$mem)[ApproxTruth[[1]]],ApproxTruth[[2]])

# run t-SNE, colored by clustering label
clust<-unlist(result$mem)
tsne_clust <- data.frame(Barcode=tsne_result[,1],TSNE.1=tsne_result[,2],
                         TSNE.2=tsne_result[,3],Clust=clust)
pdf(file="HumanSkin_BAMMSC.pdf")
visualize_clusters(as.character(tsne_clust$Clust),tsne_clust[c("TSNE.1","TSNE.2")],
                   title="BAMM-SC clustering",marker_size=0.8)
dev.off()

# run t-SNE, colored by Sample ID #
batch<-c(rep("Sample 1",dim(count[[1]])[2]),rep("Sample 2",dim(count[[2]])[2]),
         rep("Sample 3",dim(count[[3]])[2]),rep("Sample 4",dim(count[[4]])[2]),
         rep("Sample 5",dim(count[[5]])[2]))
tsne_clust <- data.frame(Barcode=tsne_result[,1],TSNE.1=tsne_result[,2],
                         TSNE.2=tsne_result[,3],Batch=batch)
pdf(file="HumanSkin_Sample.pdf")
visualize_clusters(as.character(tsne_clust$Batch),tsne_clust[c("TSNE.1","TSNE.2")],
                   title="Sample Labels",marker_size=0.8)
dev.off()

```

The t-SNE plot of clustering label is shown as:
<img src="figures/HumanSkin_BAMMSC.png" style="display: block; margin: auto;" />


The t-SNE plot of sample ID is shown as:
<img src="figures/HumanSkin_Sample.png" style="display: block; margin: auto;" />
