est.shift<-function(sample1,sample2,B=1000,min.iter=0,batch=10,mcmc.obj=NULL,dye.swap=FALSE,nb.col1=NULL, all.out=TRUE)
  {



###  Only take the finite observations
    indf1<-is.row.na(sample1)
    indf2<-is.row.na(sample2)
    sample1<-as.matrix(sample1[indf1 & indf2,])
    sample2<-as.matrix(sample2[indf1 & indf2,])
    n<-dim(sample1) 
    
    if(dye.swap==TRUE)
      {
        if(length(nb.col1)==0)
          stop("No value has been set for nb.col1" , call. = TRUE)
        else if(nb.col1>=n[2] | nb.col1<=0)
          stop("nb.col1 should be at least 0 and at most n" , call. = TRUE)
      }
    else ## just set a value for the code to run
      nb.col1<-n[2]/2 
    
    

    vec1<-as.double(t(sample1))
    vec1[is.finite(vec1)==FALSE]<- -9999999
    vec2<-as.double(t(sample2))
    vec2[is.finite(vec2)==FALSE]<- -9999999
    
    df.choice<-c(1:10,seq(20,100,10))
    df.in<-rep(100,n[2])
    w.in<-rep(1,n[1]*n[2])
        
    
    ## Minimum shift to make all the data >0
    m1<-max(0,-min(vec1)+0.01)
    m2<-max(0,-min(vec2)+0.01)
    
    
    if(length(mcmc.obj)>0)
      {
        if(class(mcmc.obj)!="mcmc.shift")
          stop("'mcmc.obj' should be of type 'mcmc.shift'" , call. = TRUE)
        
        n.iter<-length(mcmc.obj$mu)

        lambda.eps1<-mcmc.obj$lambda.eps1[n.iter]
        lambda.eps2<-mcmc.obj$lambda.eps2[n.iter]
        lambda.gamma1<-mcmc.obj$lambda.gamma1[n.iter]
        lambda.gamma2<-mcmc.obj$lambda.gamma2[n.iter]
        
        rho<-mcmc.obj$rho[n.iter]

        mu<-mcmc.obj$mu[n.iter]
        beta2<-mcmc.obj$beta[n.iter]
        alpha2<-mcmc.obj$alpha[n.iter]
        gamma1<-mcmc.obj$gamma1[,n.iter]
        gamma2<-mcmc.obj$gamma2[,n.iter]
        delta22<-mcmc.obj$delta22[n.iter]
        eta<-mcmc.obj$eta[,n.iter]
        shift<-mcmc.obj$shift[n.iter]
      }
    else
      {
        shift<-max(m1,m2)+10
        obj<-ls.effect(log2(sample1+shift),log2(sample2+shift),dye.swap=TRUE,nb.col1=nb.col1)
        mu<-obj$mu
        alpha2<-obj$alpha2
        gamma1<-obj$gamma1
        gamma2<-obj$gamma2
        if(dye.swap)
          {
            beta2<-obj$beta2
            delta22<-obj$delta22
          }
        else
          {
            beta2<-obj$beta2
            delta22<-obj$delta22
          }
        eta<-obj$eta
        if(n[2]>1)
          {
            lambda.eps1<-1/mean(obj$R1^2)
            lambda.eps2<-1/mean(obj$R2^2)
          }
        else
          {
            lambda.eps1<-1.
            lambda.eps2<-1.
          }
        
        lambda.gamma1<-0.5
        lambda.gamma2<-0.5        
        rho<-0
        
      }


    
### Main code linked to a c function
        
        
    if(all.out==TRUE)
      length<-(B-min.iter)/batch
    else
      length<-1 
   
    obj<-.C("R_link_mcmc_shift",
            vec1,
            vec2,
            as.integer(n[1]),
            as.integer(n[2]),
            as.integer(nb.col1),
            as.integer(dye.swap),
            as.integer(B),
            as.double(gamma1),
            gamma1=double(length*n[1]),
            as.double(gamma2),
            gamma2=double(length*n[1]),
            as.double(mu),
            mu=double(length),
            as.double(beta2),
            beta2=double(length),
            as.double(alpha2),
            alpha2=double(length),
            as.double(delta22),
            delta22=double(length),
            as.double(eta),
            eta=double(n[2]*length),
            as.double(lambda.eps1),
            lambda.eps1=double(length),
            as.double(lambda.eps2),
            lambda.eps2=double(length),
            as.double(w.in),
            as.double(df.choice),
            as.integer(length(df.choice)),
            as.double(df.in),
            as.double(rho),
            rho=double(length),
            as.double(shift),
            shift=double(length),
            as.double(max(-shift,0)+0.01),
            as.double(lambda.gamma1),
            as.double(lambda.gamma2),
            lambda.gamma1=double(length),
            lambda.gamma2=double(length),
            as.integer(min.iter),
            as.integer(batch),
            as.integer(all.out))

    new.mcmc<-list(gamma1=t(matrix(obj$gamma1,length,n[1],byrow=TRUE)),
                   gamma2=t(matrix(obj$gamma2,length,n[1],byrow=TRUE)),
                   mu=obj$mu, beta2=obj$beta2, alpha2=obj$alpha2,delta22=obj$delta22,
                   lambda.eps1=obj$lambda.eps1,
                   lambda.eps2=obj$lambda.eps2,
                   lambda.gamma1=obj$lambda.gamma1,lambda.gamma2=obj$lambda.gamma2,rho=obj$rho,
                   shift=obj$shift,eta=t(matrix(obj$eta,length,n[2],byrow=TRUE)))
    
### Give it the right class
    class(new.mcmc)<-"mcmc.shift"

    return(new.mcmc)
}

