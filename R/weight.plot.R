weight.plot<-function(mcmc.obj,coordinate,array=1)
  {
    l<-length(mcmc.obj$w[,array])
    
### Check the dimension matching
    if(((dim(coordinate)[1])!=l) | ((dim(coordinate)[2])!=2))
      stop("Error in the dimensions of coordinate! " , call. = TRUE)

    if(class(mcmc.obj)!="mcmc")
      stop("'mcmc.obj' should be of type 'mcmc'" , call. = TRUE)
    
    ordered.weight<-reorder.row(cbind(coordinate,mcmc.obj$w[,array]))

### Number of rows and columns
    n1<-max(coordinate[,1]+1)
    n2<-max(coordinate[,2]+1)
    
    mat.weight<-matrix(ordered.weight[,3],n1,n2,byrow=TRUE)
    title<-cbind("array", as.character(array))
    image(1:n1,1:n2,mat.weight,xlab="Row", ylab="Column",col=grey(100:0/100),main=c("Image plot of the weights ",title))
  }
