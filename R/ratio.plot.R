ratio.plot<-function(mcmc.obj,col=1,pch=1)
  {
    if(class(mcmc.obj)!="mcmc")
      stop("'mcmc.obj' should be of type 'mcmc'" , call. = TRUE)
    
    
    ### Estimate the gene effects in sample 1
    x1<-mat.mean(mcmc.obj$gamma1)[,1]
    x2<-mat.mean(mcmc.obj$gamma2)[,1]

    plot((x1+x2)/2,x1-x2,pch=pch,col=col,xlab="Overall intensity (log2)",ylab="Log ratio (log2)")
    
  }
