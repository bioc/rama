mat.mean<-function(data)
  {
    ### Main code linked to a c function
    n<-dim(data)
    vec<-as.double(t(data))
    vec[is.finite(vec)==FALSE]<- -9999999
    
    
    obj<-.C("link_R_mean_sd",
            as.double(vec),
            as.integer(n[1]),
            as.integer(n[2]),
            mean=double(n[1]),
            sd=double(n[1]))

    mean.data<-obj$mean
    sd.data<-obj$sd
    mean.data[mean.data==-9999999]<- NA
    sd.data[sd.data==-9999999]<- NA
    
    cbind(mean.data,sd.data)
  }
            
            
