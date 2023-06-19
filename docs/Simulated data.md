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

3. Apply MDS, tSNE, UMAP
```
res_MDS=cmdscale(dist(ma))
colnames(res_MDS) <- c("First Dimension", "Second Dimension")
res_tSNE=Rtsne(ma)$Y
colnames(res_tSNE) <- c("First Dimension", "Second Dimension")
res_UMAP = umap(ma)$layout
colnames(res_UMAP) <- c("First Dimension", "Second Dimension")
```
4. Apply KODAMA
```
kk=KODAMA.matrix(ma)
res_KODAMA_MDS=KODAMA.visualization(kk,method = "MDS")
res_KODAMA_tSNE=KODAMA.visualization(kk,method = "t-SNE")
res_KODAMA_UMAP=KODAMA.visualization(kk,method = "UMAP")
```

5. Plot the results

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

#### Compareing efficiency of clustering algorithm with different noisy 
```
sum <- data.frame(matrix(NA, nrow=20))
df <- data.frame()
files <- list.files(path=address, pattern="*.xlsx")
for (i in 1:length(files) {
 name <- tools::file_path_sans_ext(files[i])
 file <- as.data.frame(read_excel(files[i])) ## if you have headers in your files ##
 df <- as.data.frame(t(as.data.frame(sapply(file,CI))))
 rownames(df) <- paste(rep(1:20),name, sep = ":")
 sum <- cbind(sum, df)
 names[i] <- name }
sum <- sum[,-1]
sum <- melt(setDT(sum), measure.vars = patterns("^mean", "^lower", "^upper"),
             value.name = c("COEF_EST", "COEF_LOWER", "COEF_UPPER"))
noisy= rep(c(1:20),6)
sum <- cbind(noisy, sum)
colnames(sum)[colnames(sum)=="variable"] <- "test"
ggplot(sum, aes(x = noisy, y = COEF_EST, ymin = COEF_LOWER, ymax = COEF_UPPER)) +
  geom_ribbon(aes(fill = test), alpha = 0.3) +
  geom_line(aes(color = test))
```
<p>
  <p align="center">
    <img src="https://github.com/tkcaccia/KODAMA/blob/main/figures/CI.png" alt="hello-light" height="500" width="700" />
  </p>
</p>

