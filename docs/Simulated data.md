## Example 1: Simulated data set 

KODAMA, tSNE, UMAP are applied to simulated data set of two dimention with different degrees of noisy (from 0 to 20). The following script to compare the effect of different dimentionallity reduction algorithms on a simulated data set of 8 noisy dimensions.

### Tutorial
#### Required libraries
```
library(ggplot2)
library(gridExtra)
library(Rmisc)
library(gmodels)
```

The data is simulated with vertix function from KODAMA package with 2 dimensions and 8 noisy dimensions and then scaled

```
ma <- vertex(c(1,10), dims = 2, noisy_dimension = 8, size_cluster = 50)
```
Apply MDS, tSNE, UMAP on simulated dataset
```
res_MDS=cmdscale(dist(ma))
colnames(res_MDS) <- c("First Dimension", "Second Dimension")
res_tSNE=Rtsne(ma)$Y
colnames(res_tSNE) <- c("First Dimension", "Second Dimension")
res_UMAP = umap(ma)$layout
colnames(res_UMAP) <- c("First Dimension", "Second Dimension")
```
### KODAMA
```
kk=KODAMA.matrix(ma)
res_KODAMA_MDS=KODAMA.visualization(kk,method = "MDS")
res_KODAMA_tSNE=KODAMA.visualization(kk,method = "t-SNE")
res_KODAMA_UMAP=KODAMA.visualization(kk,method = "UMAP")
```

Visulaize the differece between different dimensionality reduction algorithm 

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

#### Generate simulated data of different noisy dimensions(1-20), then apply different algorithms and calculate the clustering efficiency of each one at different noisy levels using silhouette test.

```
noisy_dimension=c(0:20)
dimensions=2
size_cluster=50
cluster_number=2^dimensions
MDS <- matrix(nrow = 100,ncol = 20)
KODAMA_MDS <- matrix(nrow = 100,ncol = 20)
TSNE <- matrix(nrow = 100,ncol = 20)
KODAMA_TSNE <- matrix(nrow = 100,ncol = 20)
UMAP <- matrix(nrow = 100,ncol = 20)
KODAMA_UMAP <- matrix(nrow = 100,ncol = 20)

#generate simulated data
k=1
for (k in k:100){
  print(paste0 ("I am analyzing: ", k))
  for (i in 0:20 ){
    print(paste("I am anlayzing noisy", i))
    #par(mar = c(3,3,3,3))
    #par(mfrow = c(2,3))    
    ma <- vertex(c(1,10), dims = 2, noisy_dimension = i, size_cluster = 50)
    # plot(ma,pch=21,bg=rep(2:5,each=50),main="Data")
    res_MDS=cmdscale(dist(ma))
    res_tSNE=Rtsne(ma,perplexity = 20)$Y
    custom.settings = umap.defaults
    custom.settings$n_neighbors=20
    res_UMAP = umap(ma, config = custom.settings)$layout
    
    sil1 <- round(summary(silhouette(rep(1:4,each=50),dist(res_MDS)))$si.summary[4], digit=5)
    sil2 <- round(summary(silhouette(rep(1:4,each=50),dist(res_tSNE)))$si.summary[4], digit=5)
    sil3 <- round(summary(silhouette(rep(1:4,each=50),dist( res_UMAP)))$si.summary[4], digit=5)
    
    kk= KODAMA.matrix(ma)
    res_KODAMA_MDS=KODAMA.visualization(kk,method = "MDS")
    
    custom.settings = Rtsne.defaults
    custom.settings$perplexity=20
    res_KODAMA_tSNE=KODAMA.visualization(kk,method = "t-SNE",config = custom.settings)
    
    custom.settings = umap.defaults
    custom.settings$n_neighbors=20
    res_KODAMA_UMAP=KODAMA.visualization(kk,method = "UMAP",config = custom.settings)
    
    sil4 <- round(summary(silhouette(rep(1:4,each=50),dist(res_KODAMA_MDS)))$si.summary[4], digit=5)
    sil5 <- round(summary(silhouette(rep(1:4,each=50),dist(res_KODAMA_tSNE)))$si.summary[4], digit=5)
    sil6 <- round(summary(silhouette(rep(1:4,each=50),dist( res_KODAMA_UMAP)))$si.summary[4], digit=5)
    
    mtext(paste0("Noisy",i), side = 3, line = -25, outer = TRUE)
    
    MDS [k,i] <- sil1
    TSNE [k,i] <- sil2
    UMAP [k,i] <- sil3
    
    KODAMA_MDS [k,i] <- sil4
    KODAMA_TSNE [k,i]<- sil5
    KODAMA_UMAP [k,i]<- sil6
  }
}
```
####  The confidence intervals for each clustering algorithm at different noisy level  are calculated and visualized 

```
test <-  list(MDS, TSNE,UMAP,K_MDS,K_TSNE,K_UMAP)
names(test) <- c("MDS", "TSNE","UMAP","KODAMA_MDS","KODAMA_TSNE","KODAMA_UMAP")
sum <- data.frame(matrix(NA, nrow=20))
df <- data.frame()
for (i in 1:length(test)) {
  df <- as.data.frame(t(as.data.frame(sapply(test[[i]],CI))))
  sum <- cbind(sum, df) }
  sum <- sum[,-1]
  sum <- melt(setDT(sum), measure.vars = patterns("^mean", "^lower", "^upper"),
              value.name = c("COEF_EST", "COEF_LOWER", "COEF_UPPER"))
  noisy= rep(c(1:20),6)
  sum <- cbind(noisy, sum)
  sum$test <- rep(names(test), each= 20)
  
  ggplot(sum, aes(x = noisy, y = COEF_EST, ymin = COEF_LOWER, ymax = COEF_UPPER)) +
    geom_ribbon(aes(fill = test), alpha = 0.3) +
    geom_line(aes(color = test))+
    ylab("95% confidence inteval of silhouette scores")+
    xlab("Noisy dimensions")
```
<p>
  <p align="center">
    <img src="https://github.com/tkcaccia/KODAMA/blob/main/figures/CI%20simulated.png" alt="hello-light" />
  </p>
</p>

