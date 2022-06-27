# KODAMA

This is a version in developing of KODAMA with landmarks to adapt the algorithm to the analysis of data set with more than 10,000 entries.
The wrapper for the C++ implementation of Barnes-Hut t-Distributed Stochastic Neighbor Embedding has been integrated to convert the KODAMA's dissimilarity matrix in a low dimensional space.

## Metabolomic data

The data belong to a cohort of 22 healthy donors (11 male and 11 female) where each provided about 40 urine samples over the time course of approximately 2 months, for a total of 873 samples. Each sample was analysed by Nuclear Magnetic Resonance Spectroscopy. Each spectrum was divided in 450 spectral bins.

```
data(MetRef)
u=MetRef$data;
u=u[,-which(colSums(u)==0)]
u=normalization(u)$newXtrain
u=scaling(u)$newXtrain
class=as.numeric(as.factor(MetRef$gender))
plot(cc,pch=21,bg=class,xlab="First Component",ylab="Second Component")

class=as.numeric(as.factor(MetRef$donor))
kk=KODAMA(u,landmarks = 100,perplexity = 10,f.par = 50)
tt=Rtsne(u,perplexity = 10)

par(mfrow=c(1,2))
plot(kk$scores,pch=21,bg=rainbow(22)[class], main="KODAMA-tSNE")
plot(tt$Y,pch=21,bg=rainbow(22)[class], main="tSNE, xlab="First Component", ylab="Second Component")
```
![This is an image](https://github.com/tkcaccia/Documents/blob/main/MetRef.png)


## Single-cell data

```


ta=read.csv("tasic-sample_heatmap_plot_data.txt")
rownames(ta)=ta[,1]
VIS=read.csv("mouse_VISp_gene_expression_matrices_2018-06-14/mouse_VISp_2018-06-14_exon-matrix.csv")
ALM=read.csv("mouse_ALM_gene_expression_matrices_2018-06-14/mouse_ALM_2018-06-14_exon-matrix.csv")

data=t(cbind(ALM,VIS))
colnames(data)=as.character(data[1,])

data=data[-1,]

ii=intersect(rownames(data),rownames(ta))
data=data[ii,]

data=data[,colSums(data)!=0]

near.zero.counts=colMeans(data<32)
temp=data
temp[temp<=32]=NA
temp=log2(temp)
m=colMeans(temp,na.rm = TRUE)

y=exp(-1.5*(m-6.56))+0.02


data=data[,which(near.zero.counts>y)]

su=rowSums(data)
data=((data/su)*10^6)*median(su)
data=log2(data+1)

pca=prcomp(data)$x[,1:50]

kk=KODAMA(pca,landmarks = 1000)
plot(kk$scores,pch=21,bg=ta[,"cluster_color"])

tt=Rtsne(pca)
```


# GeoMX

https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.R
https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html#3_Study_Design
