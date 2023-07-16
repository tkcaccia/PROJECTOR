set.seed(1)
n <- 10
x <- runif(n)
y <- runif(n)
p <- cbind(x,y)
xlim <- c(min(x) - 0.1*diff(range(x)), c(max(x) + 0.1*diff(range(x))))
ylim <- c(min(y) - 0.1*diff(range(y)), c(max(y) + 0.1*diff(range(y))))
plot(p, xlim=xlim, ylim=ylim)
text(p, labels=seq(n), pos=3)


xspline(x, y, shape = c(0,rep(-1, 10-2),0), border="red")

dd= xspline(x, y, shape = c(0,rep(-1, 10-2),0), border="red",draw = FALSE)

plot(x,y)
new_trajectory = function(x,y,n=20,data=NULL,knn=10,FUN=mean){
  ii=identify(x,y,order = TRUE)
  ii=ii$ind[order(ii$order)]
  dd= xspline(x[ii], y[ii], shape = c(0,rep(-1, 10-2),0), border="red",draw = FALSE)
  ll=length(dd$x)
  sel=seq(1,ll,length.out =n)
  
 # points(dd,col=2,bg="#eeeeee",lwd=2,pch=21)
  dd$x=dd$x[sel]
  dd$y=dd$y[sel]
  xy=cbind(dd$x,dd$y)
  xy_total=cbind(x,y)
  selection=knn_Armadillo(xy_total,xy,k = knn)$nn_index
  if(!is.null(data)){
     trajectory=apply(selection,1,function(z) apply(data[z,],2,FUN))
  }
  points(dd,col=2,bg="#eeeeee",lwd=2,pch=21)
  list(xy=dd,selection=selection,trajectory=trajectory,
       settings=list(x=x,y=y,n=n,data=data,knn=knn,FUN=FUN))
}

add_branch = function(dd){
  n_start=identify(dd$xy,n=1)
  start_x=dd$xy$x[n_start]
  start_y=dd$xy$y[n_start]
  ii=identify(x,y,order = TRUE)
  ii=ii$ind[order(ii$order)]
  
  branch= xspline(c(start_x,x[ii]), c(start_y,y[ii]), shape = c(0,rep(-1, 10-2),0), border="red",draw = FALSE)
  ll=length(branch$x)
  sel=seq(1,ll,length.out =dd$settings$n-n_start+1)
  
  # points(dd,col=2,bg="#eeeeee",lwd=2,pch=21)
  branch$x=branch$x[sel]
  branch$y=branch$y[sel]
  
  
  xy=cbind(branch$x,branch$y)
  xy_total=cbind(dd$settings$x,dd$settings$y)
  if(!is.null(data)){
    selection=knn_Armadillo(xy_total,xy,k = dd$settings$knn)$nn_index
    trajectory=apply(selection,1,function(z) apply(data[z,],2,FUN))
  }
  points(branch,col=3,bg="#eeeeee",lwd=2,pch=21)
  
}

dd=new_trajectory(x,y,20)
add_branch(dd)
