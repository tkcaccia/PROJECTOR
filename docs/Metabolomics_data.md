## Metabolomic data

The data belong to a cohort of 22 healthy donors (11 male and 11 female) where each provided about 40 urine samples over the time course of approximately 2 months, for a total of 873 samples. Each sample was analysed by Nuclear Magnetic Resonance Spectroscopy. Each spectrum was divided in 450 spectral bins.

### Tutorial

1. Data upload and processing 

```
data(MetRef)
u=MetRef$data
u=u[,-which(colSums(u)==0)]
u=normalization(u)$newXtrain
u=scaling(u)$newXtrain
class=as.numeric(as.factor(MetRef$gender))
class2=as.numeric(as.factor(MetRef$donor))
```

2. Apply MDS, tSNE and UMAP
```
res_MDS=cmdscale(dist(u))
res_tSNE=Rtsne(u)$Y
res_UMAP = umap(u)$layout
```

3. Apply KODAMA (The f. par has been adjusted to 50 which means 50 princible components for the Partial Least Square "PLS" classifier)

```
kk=KODAMA.matrix(u,f.par = 50)
res_KODAMA_MDS=KODAMA.visualization(kk,method = "MDS")
res_KODAMA_tSNE=KODAMA.visualization(kk,method = "t-SNE")
res_KODAMA_UMAP=KODAMA.visualization(kk,method = "UMAP")
```

4. Visualize the different clustering algorithmss:

  a) labelled by the gender

```
par(mfrow = c(2,3))
plot(res_MDS,pch=21,bg=rainbow(2)[class],main="MDS")
plot(res_tSNE,pch=21,bg=rainbow(2)[class],main="tSNE")
plot(res_UMAP,pch=21,bg=rainbow(2)[class],main="UMAP")
plot(res_KODAMA_MDS,pch=21,bg=rainbow(2)[class],main="KODAMA_MDS",)
plot(res_KODAMA_tSNE,pch=21,bg=rainbow(2)[class],main="KODAMA_tSNE")
plot(res_KODAMA_UMAP,pch=21,bg=rainbow(2)[class],main="KODAMA_UMAP")
```
<p>
  <p align="center">
    <img src="https://github.com/ebtesam-rashid/KODAMA.Caccio/blob/main/Figures/metab%20gender%20final.png" alt="hello-light" height="500" width="700" />
  </p>
</p>

  b) labelled by the donor

```
plot(res_MDS,pch=21,bg=rainbow(22)[class2],main="MDS")
plot(res_tSNE,pch=21,bg=rainbow(22)[class2],main="tSNE")
plot(res_UMAP,pch=21,bg=rainbow(22)[class2],main="UMAP")
plot(res_KODAMA_MDS,pch=21,bg=rainbow(22)[class2],main="KODAMA_MDS",)
plot(res_KODAMA_tSNE,pch=21,bg=rainbow(22)[class2],main="KODAMA_tSNE")
plot(res_KODAMA_UMAP,pch=21,bg=rainbow(22)[class2],main="KODAMA_UMAP")
```
<p>
  <p align="center">
    <img src="https://github.com/ebtesam-rashid/KODAMA.Caccio/blob/main/Figures/final%20metab%20donor.png" alt="hello-light" height="500" width="700" />
  </p>
</p>

