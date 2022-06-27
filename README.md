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
