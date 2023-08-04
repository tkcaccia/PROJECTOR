## Simulated data sets

In this code, simulated data sets with a different number of noisy dimensions ranging from 0 to 20 are generated. For each number of noisy dimensions, 100 different data sets are generated. The results of each method are compared using the Silhouette function.

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
  print(paste0 ("Repetition number: ", k))
  for (i in 0:20 ){
    print(paste("Number of noisy dimensions:", i))
    ma <- vertex(c(1,10), dims = 2, noisy_dimension = i, size_cluster = 50)
    res_MDS=cmdscale(dist(ma))
    res_tSNE=Rtsne(ma)$Y
    res_UMAP = umap(ma)$layout
    
    sil1 <- round(summary(silhouette(rep(1:4,each=50),dist(res_MDS)))$si.summary[4], digit=5)
    sil2 <- round(summary(silhouette(rep(1:4,each=50),dist(res_tSNE)))$si.summary[4], digit=5)
    sil3 <- round(summary(silhouette(rep(1:4,each=50),dist( res_UMAP)))$si.summary[4], digit=5)
    
    
    kk=KODAMA.matrix(ma,FUN="KNNPLS-DA",spatial.knn = 10)
    res_KODAMA_MDS=KODAMA.visualization(kk,method = "MDS")
    res_KODAMA_tSNE=KODAMA.visualization(kk,method = "t-SNE")
    res_KODAMA_UMAP=KODAMA.visualization(kk,method = "UMAP")
    
    sil4 <- round(summary(silhouette(rep(1:4,each=50),dist(res_KODAMA_MDS)))$si.summary[4], digit=5)
    sil5 <- round(summary(silhouette(rep(1:4,each=50),dist(res_KODAMA_tSNE)))$si.summary[4], digit=5)
    sil6 <- round(summary(silhouette(rep(1:4,each=50),dist( res_KODAMA_UMAP)))$si.summary[4], digit=5)
    

    MDS [k,i+1] <- sil1
    TSNE [k,i+1] <- sil2
    UMAP [k,i+1] <- sil3
    
    KODAMA_MDS [k,i+1] <- sil4
    KODAMA_TSNE [k,i+1]<- sil5
    KODAMA_UMAP [k,i+1]<- sil6
  }
}

```
#### Required libraries
```
library(ggplot2)
library(cluster)
library(gridExtra)
library(Rmisc)
library(gmodels)
```
The confidence intervals for each clustering algorithm at different noisy level  are calculated and visualized 

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

