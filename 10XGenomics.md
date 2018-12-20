# BAMMSC for 10X Genomics

The output of 10X Genomics cellrangerRkit consists of three files that are matrix.mtx, genes.tsv and barcodes.tsv. In this
tutorial, we demonstrate how to process the three files to obtain the input format that is needed by BAMMSC for clustering.

**Step1: Create a data matrix of for each invidual from output of 10X Genomics cellrangerRkit**
----------------------
```
library(Matrix)
mat=readMM("matrix.mtx") 
gene_info=read.delim("genes.tsv", stringsAsFactors=FALSE, sep="\t", header=FALSE) 
barcodes=read.delim("barcodes.tsv", stringsAsFactors=FALSE, sep="\t", header=FALSE) 
rownames(mat)= gene_info[, 1] 
colnames(mat)= barcodes[,1] 
mat = as.matrix(mat) 
pd = data.frame(id = barcodes[, 1], row.names = barcodes[, 1]) 
colnames(pd) = c("barcode") 
gene_symbols = gene_info 
row.names(gene_symbols) = gene_info[, 1] 
colnames(gene_symbols) = c("id", "symbol") 
gbm= newGeneBCMatrix(mat = mat, fd = gene_symbols, pd = pd) 
gbm@barcode_filtered= TRUE

#For individual 1, create raw count matrix
data_sample_1=as.matrix(exprs(gbm))
rownames(data_sample_1)=fData(gbm)$symbol

# Repeat the above steps for each individual to get data_sample_1, data_sample_2, ... , data_sample_L
```


**Step2: Merge mutiple data matrix together to form a list of data matrix**
----------------------
### Create a list to include all individuals and generate input for BAMM-SC
```
count=list()
count[[1]]=data_sample_1
count[[2]]=data_sample_2
count[[3]]=data_sample_3
          .
          .
          .
count[[L]]=data_sample_L
```

**Step3: Gene filtering by variance before running BAMMSC**
----------------------
```
countall=do.call(cbind,count)
genesd=apply(countall,1,sd)
ngene=1000
select=sort(order(genesd)[dim(countall)[1]:1][1:ngene])
for(l in 1:length(count)){
  count[[l]]=count[[l]][select,]
}
```


**Step4: Run BAMMSC for cell clusetering**
----------------------
```
library(BAMMSC)
result=BAMMSC(data, K=4)
```
