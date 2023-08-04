# KODAMA for Spatial Transcriptomics

We compared the output of KODAMA with other dimensionality reduction methods such as Uniform Manifold Approximation and Projection (UMAP) and t-Distributed Stochastic Neighbour Embedding (t-SNE).

## Example 1: Simulated data set 
A simulated dataset is generated distributing the values of two dimensions centered in the vertix of a square. Additional non-informative dimensions (noise) are included. In the following script, we compare the results of different dimensionality reduction algorithms on a simulated data set with 8 noisy dimensions.

The data is simulated with vertex function from KODAMA package with 2 dimensions and 8 noisy dimensions and then scaled

```
ma <- vertex(c(1,10), dims = 2, noisy_dimension = 8, size_cluster = 50)
```

We perform the analysis of MDS, t-SNE, and UMAP.
```
res_MDS=cmdscale(dist(ma))
colnames(res_MDS) <- c("First Dimension", "Second Dimension")
res_tSNE=Rtsne(ma)$Y
colnames(res_tSNE) <- c("First Dimension", "Second Dimension")
res_UMAP = umap(ma)$layout
colnames(res_UMAP) <- c("First Dimension", "Second Dimension")
```
The KODAMA analysis is performed and the KODAMA's dissimilarity matrix is transformed in a two-dimensional space using MDS, t-SNE, and UMAP. 
```
kk=KODAMA.matrix(ma,FUN="KNNPLS-DA",spatial.knn = 10)
res_KODAMA_MDS=KODAMA.visualization(kk,method = "MDS")
res_KODAMA_tSNE=KODAMA.visualization(kk,method = "t-SNE")
res_KODAMA_UMAP=KODAMA.visualization(kk,method = "UMAP")
```

The results are plotted.

```
par(mfrow = c(2,3))
labels <- rep(c("#FF0000","#0000FF","#008000","#FFFF00"),each= 50)
plot(res_MDS,pch=21,bg=labels,main="MDS")
plot(res_tSNE,pch=21,bg=labels,main="tSNE")
plot(res_UMAP,pch=21,bg=labels,main="UMAP")
plot(res_KODAMA_MDS,pch=21,bg=labels,main="KODAMA_MDS",ylim=range(res_KODAMA_MDS[,1]))
plot(res_KODAMA_tSNE,pch=21,bg=labels,main="KODAMA_tSNE")
plot(res_KODAMA_UMAP,pch=21,bg=labels,main="KODAMA_UMAP")


```
<p>
  <p align="center">
    <img src="https://github.com/tkcaccia/KODAMA/blob/main/figures/Rplot02.png" alt="hello-light" />
  </p>
</p>

We compared now the output Simulated data of different noisy dimensions(1-20) are generated. Then apply different algorithms and calculate the clustering efficiency of each one at different noisy levels using silhouette test.The confidence intervals for each clustering algorithm at different noisy level  are calculated and visualized. [Simulated data](https://github.com/tkcaccia/KODAMA/edit/main/docs/Simulated%20data.md). The clustering quality of KODAM is high compared to othe ther algorithms. 

<p>
  <p align="center">
    <img src="https://github.com/tkcaccia/KODAMA/blob/main/figures/CI%20simulated.png" alt="hello-light" height="500" width="700" />
  </p>
</p>

## Example 2: GEOMx dataset 1
The GeoMx Digital Spatial Profiler (DSP) is a platform for capturing spatially resolved high-plex gene (or protein) expression data from tissue [Merritt et al., 2020](https://pubmed.ncbi.nlm.nih.gov/32393914/). In particular, formalin-fixed paraffin-embedded (FFPE) or fresh-frozen (FF) tissue sections are stained with barcoded in-situ hybridization probes that bind to endogenous mRNA transcripts. 
GeoMx kidney dataset has been created with the human whole transcriptome atlas (WTA) assay. The dataset includes 4 diabetic kidney disease (DKD) and 3 healthy kidney tissue samples. Regions of interest (ROI) were spatially profiled to focus on two different kidney structures: tubules or glomeruli. One glomerular ROI contains the entirety of a single glomerulus. Each tubular ROI contains multiple tubules that were segmented into distal (PanCK+) and proximal (PanCK-) tubule areas of illumination (AOI). The preprocessing workflow is described [here](https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html).
An imputing procedure was added to the original [R script](https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.R).

### Tutorial

#### Install required packages
```
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("impute")
install.packages("KODAMA")
library(impute)
library(KODAMA)
```
#### Ulpoad data
```
data=t(log2(assayDataElement(target_demoData , elt = "q_norm")))
data[is.infinite(data)]=NA
data=impute.knn(data)$data
```
#### Run MDS
```
MDS_out=cmdscale(dist(data))
pData(target_demoData)[, c("MDS1", "MDS2")] <- MDS_out[, c(1,2)]
```
#### run tSNE
```
set.seed(42) # set the seed for tSNE as well
tsne_out <- Rtsne(data, perplexity = ncol(target_demoData)*.15)
pData(target_demoData)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
```
#### run UMAP
```
umap_out <- umap(data, config = custom_umap)
pData(target_demoData)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
```
#### run KODAMA
```
kk=KODAMA.matrix(data)
res= KODAMA.visualization(kk)
res1= KODAMA.visualization(kk,method = "MDS")
res2= KODAMA.visualization(kk,method = "t-SNE")
res3= KODAMA.visualization(kk,method = "UMAP")
pData(target_demoData)[, c("KODAMA1.MDS", "KODAMA2.MDS")] <- res1
pData(target_demoData)[, c("KODAMA1.tSNE", "KODAMA2.tSNE")] <- res2
pData(target_demoData)[, c("KODAMA1.UMAP", "KODAMA2.UMAP")] <- res3
```
#### MDS vs KODAMA.MDS
```
plot1=ggplot(pData(target_demoData), aes(x = MDS1, y = MDS2,color = segment, shape = class)) + geom_point(size = 3) + theme_bw()
plot2=ggplot(pData(target_demoData), aes(x = KODAMA1.MDS, y = KODAMA2.MDS, color = segment, shape = class)) + geom_point(size = 3) + theme_bw()
grid.arrange(plot1, plot2, ncol=2)
```
<p>
  <p align="center">
    <img src="https://github.com/tkcaccia/KODAMA/blob/main/figures/MDS%20geomx.png" alt="hello-light" height="500" width="700" />
  </p>
</p>

#### tSNA vs KODAMA.tSNE
```
plot3=ggplot(pData(target_demoData), aes(x = tSNE1, y = tSNE2, color = segment, shape = class)) + geom_point(size = 3) + theme_bw()
plot4=ggplot(pData(target_demoData), aes(x = KODAMA1.tSNE, y = KODAMA2.tSNE, color = segment, shape = class)) + geom_point(size = 3) + theme_bw()
grid.arrange(plot3, plot4, ncol=2)
```
<p>
  <p align="center">
    <img src="https://github.com/tkcaccia/KODAMA/blob/main/figures/tsne%20geomx.png" alt="hello-light" height="500" width="700" />
  </p>
</p>

#### UMAP vs KODAMA.UMAP
```
plot5=ggplot(pData(target_demoData), aes(x = UMAP1, y = UMAP2, color = segment, shape = class)) + geom_point(size = 3) + theme_bw()
plot6=ggplot(pData(target_demoData), aes(x = KODAMA1.UMAP, y = KODAMA2.UMAP, color = segment, shape = class)) + geom_point(size = 3) + theme_bw()
grid.arrange(plot5, plot6, ncol=2)
```
<p>
  <p align="center">
    <img src="https://github.com/tkcaccia/KODAMA/blob/main/figures/umap%20geomx.png" alt="hello-light" height="500" width=700" />
  </p>
</p>

## Example 3: GeoMx dataset2
This dataset represent the quantitate transcript and protein abundance in spatially distinct regions of metastic prostate cancer tissue retreived by GeoMx Digital Spatial Profiler (DSP)  [Brady et al., 2021](https://www.nature.com/articles/s41467-021-21615-4). This study entails leaving 141 ROIs from 53 metastases (26 patients). 
                                                                                                                                           
### Tutorial
#### Install required packages
 ```
library(KODAMA)
library(readxl)
library(tidyr)
library(ggplot2)
library(gridExtra)
```
#### Ulpoad data
```
dat3 <- as.data.frame(read.csv("Supplementary_Data_File_3.txt", header=TRUE, sep = "\t", dec = "."))
dat4 <- as.data.frame(read.csv("Supplementary_Data_File_4.txt", header=TRUE, sep = "\t", dec = "."))
dat5<- as.data.frame(read_excel("41467_2021_21615_MOESM5_ESM.xlsx", skip = 1))
```
### Data preprocessing

```
rownames(dat5) <- dat5[,"Sample_ID"]
sel <- c( "Gene", "Negative_Normalized", "Sample_ID")
g <- dat3[,sel]
g2 <- as.data.frame(g |> pivot_wider(names_from = Gene, values_from =Negative_Normalized ))
rownames(g2) <- g2$Sample_ID
g2 <- g2[,-1]
sel=intersect(rownames(g2),rownames(dat5))
g3=g2[sel,]
metadata=dat5[sel,]
sel <- c( "protein", "count_ngs_norm", "Sample_ID")
Sample_ID=paste(dat4$tissue_id,dat4$punch,sep="_")
p=data.frame(Sample_ID=Sample_ID,count_ngs_norm=dat4$count_ngs_norm,Protein=dat4$Protein)
p=p[!is.na(dat4$tissue_id) & !is.na(dat4$punch),]
p2 <- as.data.frame(p |> pivot_wider(names_from = Protein, values_from =count_ngs_norm ))
rownames(p2) <- p2$Sample_ID
sel=intersect(rownames(g2),rownames(dat5))
p3=p2[sel,]
``` 
#### Run MDS
```
res_MDS=cmdscale(dist(g3))
metadata[, c("MDS1", "MDS2")] <- res_MDS[, c(1,2)]
```
#### run tSNE
```
set.seed(42) # set the seed for tSNE as well
tsne_out <- Rtsne(g3, perplexity = 10)
metadata[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
```
#### run UMAP
```
custom.settings = umap.defaults
custom.settings$n_neighbors=10
umap_out <- umap(g3,config = custom.settings)
metadata[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
```
#### run KODAMA
```
kk=KODAMA.matrix(g3)
res= KODAMA.visualization(kk)
res1= KODAMA.visualization(kk,method = "MDS")
custom.settings$perplexity=10
custom.settings = Rtsne.defaults
res2= KODAMA.visualization(kk,method = "t-SNE", config = custom.settings)
custom.settings = umap.defaults
custom.settings$n_neighbors=10
res3= KODAMA.visualization(kk,method = "UMAP", config = custom.settings)
metadata[, c("KODAMA1.MDS", "KODAMA2.MDS")] <- res1
metadata[, c("KODAMA1.tSNE", "KODAMA2.tSNE")] <- res2
metadata[, c("KODAMA1.UMAP", "KODAMA2.UMAP")] <- res3
```
#### MDS vs KODAMA.MDS
```
Histology=as.factor(metadata$`Histology category`)
plot1=ggplot(metadata, aes(x = MDS1, y = MDS2, color = class)) + geom_point(aes(fill=Histology), colour="black",pch=21, size=4) + theme_bw()
plot2=ggplot(pData(target_demoData), aes(x = KODAMA1.MDS, y = KODAMA2.MDS, color = class)) + geom_point(aes(fill=Histology), colour="black",pch=21, size=4) + theme_bw()
grid.arrange(plot1, plot2, ncol=2)
```
<p>
  <p align="center">
    <img src="https://github.com/tkcaccia/KODAMA/blob/main/figures/MDS%20dat2.png" alt="hello-light" height="500" width="700" />
  </p>
</p>

#### tSNA vs KODAMA.tSNE
```
plot3=ggplot(metadata, aes(x = tSNE1, y = tSNE2, color = class)) + geom_point(aes(fill=Histology), colour="black",pch=21, size=4) + theme_bw()
plot4=ggplot(metadata, aes(x = KODAMA1.tSNE, color = class)) + geom_point(aes(fill=Histology), colour="black",pch=21, size=4) + theme_bw()
grid.arrange(plot3, plot4, ncol=2)
```
<p>
  <p align="center">
    <img src="https://github.com/tkcaccia/KODAMA/blob/main/figures/TSNE%20dat2.png" alt="hello-light" height="500" width="700" />
  </p>
</p>

#### UMAP vs KODAM.UMAP
```
plot5=ggplot(metadata, aes(x = UMAP1, y = UMAP2, color = class)) + geom_point(aes(fill=Histology), colour="black",pch=21, size=4) + theme_bw()
plot6=ggplot(metadata, aes(x = KODAMA1.UMAP, color = class)) + geom_point(aes(fill=Histology), colour="black",pch=21, size=4) + theme_bw()
grid.arrange(plot5, plot6, ncol=2)
```
<p>
  <p align="center">
    <img src="https://github.com/tkcaccia/KODAMA/blob/main/figures/UMAP%20dat2.png" alt="hello-light" height="500" width=700" />
  </p>
</p>                  

#### KODAMA.umap clustering according to Androgen receptor(AR)

```
values=as.numeric(p3$AR)
v=quantile(values,probs=c(0.2,0.4,0.6,0.8),na.rm = TRUE)
AR.protein=findInterval(values, v)
plot7=ggplot(metadata,
             aes(x = KODAMA1umap, y = KODAMA2umap)) +
  geom_point(aes(fill=AR.protein), 
             colour="black",pch=21, size=4) +
  theme_bw()
grid.arrange(plot7, ncol=1)
```

<p>
  <p align="center">
    <img src="https://github.com/tkcaccia/KODAMA/blob/main/figures/AR.png" alt="hello-light" height="500" width=700" />
  </p>
</p> 
                                                                                                                    
#### KODAMA.umap clustering according to CD68

```                                                                                                                                                                                 
values=as.numeric(p3$CD68)
v=quantile(values,probs=c(0.2,0.4,0.6,0.8),na.rm = TRUE)
CD68=findInterval(values, v)
plot7=ggplot(metadata,
             aes(x = KODAMA1umap, y = KODAMA2umap)) +
  geom_point(aes(fill=CD68), 
             colour="black",pch=21, size=4) +
  theme_bw()
grid.arrange(plot7, ncol=1)
```
                           
<p>
  <p align="center">
    <img src="https://github.com/tkcaccia/KODAMA/blob/main/figures/CD68.png" alt="hello-light" height="400" width=700" />
  </p>
</p>                                                                                                                                 
