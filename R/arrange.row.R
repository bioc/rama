arrange.row<-function(data)
  {
### convert the data.frame to a matrix
    data<-as.matrix(data)
    n<-dim(data)

### Check that the indices start at 0
    if(min(data[,1])!=0 | min(data[,2])!=0)
      stop("The indices should start at zero" , call. = TRUE)
    
    m1<-max(data[,1]+1)
    m2<-max(data[,2]+1)
    
    obj<-.C("reorder",
            data=as.double(t(data)),
            n[1],
            n[2],
            all.data=as.double(rep(-999999,m1*m2*n[2])),
            as.integer(m1),
            as.integer(m2),PACKAGE="rama")
    
    data<-obj$all.data
    data<-data[data!=-999999]
    
    matrix(data,n[1],n[2],byrow=TRUE)
  }
