## Single-cell data

The data set from Tasic et al. encompasses 23,822 cells from adult mouse cortex, split by the authors into 133 clusters with strong hierarchical organisation. A standard preprocessing pipeline consisting of sequencing depth normalisation, feature selection, log-transformation, and reducing the dimensionality to 50 PCs was applied as described by Kobak & Berens in [The art of using t-SNE for single-cell transcriptomics](https://www.nature.com/articles/s41467-019-13056-x).

Download the data from [here](http://celltypes.brain-map.org/rnaseq) and unpack. Direct links: [VISp](http://celltypes.brain-map.org/api/v2/well_known_file_download/694413985), [ALM](http://celltypes.brain-map.org/api/v2/well_known_file_download/694413179).
To get the information about cluster colors and labels (sample_heatmap_plot_data.csv), open the interactive [data browser](http://celltypes.brain-map.org/rnaseq/mouse/v1-alm), go to "Sample Heatmaps", click "Build Plot!" and then "Download data as CSV".


### Tutorial

1. Data upload

```
ta=read.csv("tasic-sample_heatmap_plot_data.txt")
rownames(ta)=ta[,1]
VIS=read.csv("mouse_VISp_gene_expression_matrices_2018-06-14/mouse_VISp_2018-06-14_exon-matrix.csv")
ALM=read.csv("mouse_ALM_gene_expression_matrices_2018-06-14/mouse_ALM_2018-06-14_exon-matrix.csv")
```

2. Data processing

```
data=t(cbind(ALM,VIS))
colnames(data)=as.character(data[1,])
data=data[-1,]
ii=intersect(rownames(data),rownames(ta))
data=data[ii,]
data=data[,colSums(data)!=0]
near.zero.counts=colMeans(data<32)
```

3. Data normalization

```
temp=data
temp[temp<=32]=NA
temp=log2(temp)
m=colMeans(temp,na.rm = TRUE)
y=exp(-1.5*(m-6.56))+0.02
data=data[,which(near.zero.counts>y)]
su=rowSums(data)
data=((data/su)*10^6)*median(su)
data=log2(data+1)
```

4. Calculate PCA for the processed data

```
pca=prcomp(data)$x[,1:50]
```

5. Apply KODAMA

```
kk=KODAMA.matrix(pca)
res_KODAMA_tSNE <- KODAMA.visualization(kk)
plot(res_KODAMA_tSNE,pch=21,bg=ta[,"cluster_color"],main="KODAMA", xlab= "Fisrt dimension", ylab = "Second dimension")
```

6. Plot KODAMA result


<p>
  <p align="center">
    <img src="https://github.com/ebtesam-rashid/KODAMA.Caccio/blob/main/Figures/final%20single%20cell%20kodama.png" alt="hello-light" height="500" width="500" />
  </p>
</p>


