fit.model<-function(sample1,sample2,B=1000,min.iter=0,batch=10,shift=NULL,mcmc.obj=NULL, dye.swap=FALSE, nb.col1=NULL, all.out=TRUE)
  {


    ###  Only take the finite observations
    indf1<-is.row.na(sample1)
    indf2<-is.row.na(sample2)
    sample1<-as.matrix(sample1[indf1 & indf2,])
    sample2<-as.matrix(sample2[indf1 & indf2,])
    n<-dim(sample1) 
    
    vec1<-as.double(t(sample1))
    vec1[is.finite(vec1)==FALSE]<- -9999999
    vec2<-as.double(t(sample2))
    vec2[is.finite(vec2)==FALSE]<- -9999999
    
    df.choice<-c(1:10,seq(20,100,10))
    w.out<-rep(0,n[1]*n[2])
               
    a.gamma<-1
    b.gamma<-0.005
            
    if(dye.swap==TRUE)
      {
        if(length(nb.col1)==0)
          stop("No value has been set for nb.col1" , call. = TRUE)
        else if(nb.col1>=(n[2]-1) | nb.col1<=1)
          stop("nb.col1 should be at least 2 and at most n-2" , call. = TRUE)
      }
    else ## just set a value for the code to run
      nb.col1<-n[2]/2 
    
    if(length(shift)==0) ## No shift has been specified
      {
        ## Estimate the shift
        
        m1<-max(0,-min(vec1)+0.01)
        m2<-max(0,-min(vec2)+0.01)
        shift<-max(m1,m2)
        mcmc.obj.shift<-est.shift(sample1,sample2,B=50000,min.iter=10000,batch=40,mcmc.obj=NULL,dye.swap=dye.swap,nb.col1=nb.col1)
        shift<-mean(mcmc.obj.shift$shift)
      }
    
    if(length(mcmc.obj)>0)
      {
        if(class(mcmc.obj)!="mcmc")
          stop("'mcmc.obj' should be of type 'mcmc'" , call. = TRUE)
    
        
        n.iter<-length(mcmc.obj$mu)
        print(n.iter)
        
        lambda.eps1<-mcmc.obj$lambda.eps1[,n.iter]
        lambda.eps2<-mcmc.obj$lambda.eps2[,n.iter]
        lambda.gamma1<-mcmc.obj$lambda.gamma1[n.iter]
        lambda.gamma2<-mcmc.obj$lambda.gamma2[n.iter]
        
        a.eps<-mcmc.obj$a.eps[n.iter]
        b.eps<-mcmc.obj$b.eps[n.iter]
        
        rho<-mcmc.obj$rho[n.iter]

        df.in<-mcmc.obj$df[,n.iter]
        w.in<-matrix(mcmc.obj$w)
        w.in<-as.double(t(w.in))
        mu<-mcmc.obj$mu[n.iter]
        

        alpha2<-mcmc.obj$alpha[n.iter]
        gamma1<-mcmc.obj$gamma1[,n.iter]
        gamma2<-mcmc.obj$gamma2[,n.iter]
        beta2<-mcmc.obj$beta[n.iter]
        delta22<-mcmc.obj$delta22[n.iter]
        eta<-mcmc.obj$eta[,n.iter]
        
      }
    else
      {

        obj<-ls.effect(log2(sample1+shift),log2(sample2+shift),dye.swap=TRUE,nb.col1=nb.col1)
        mu<-obj$mu
        alpha2<-obj$alpha2
        
        gamma1<-obj$gamma1
        gamma2<-obj$gamma2
        if(dye.swap==TRUE)
          {
            beta2<-obj$beta2
            delta22<-obj$delta22
          }
        else
          {
            beta2<-0
            delta22<-0
          }
        eta<-obj$eta
        if(n[2]>1)
          {
            lambda.eps1<-1/mat.mean(obj$R1^2)[,1]
            lambda.eps2<-1/mat.mean(obj$R2^2)[,1]
          }
        else
          {
            lambda.eps1<-rep(1,n[1])
            lambda.eps2<-rep(1,n[1])
          }
        
        
        lambda.gamma1<-0.5
        lambda.gamma2<-0.5
        
        
        a.eps<-median(lambda.eps1)
        b.eps<-mad(lambda.eps1)
        df.in<-rep(10,n[2])
        w.in<-rep(1,n[1]*n[2])
        
        rho<-.8
                

      }

### Main code linked to a c function

    if(all.out==TRUE)
      length<-(B-min.iter)/batch
    else
      length<-1
        
    obj<-.C("ex_R_link_mcmc",
            log2(vec1+shift),
            log2(vec2+shift),
            as.integer(n[1]),
            as.integer(n[2]),
            as.integer(nb.col1),
            as.integer(B),
            as.integer(dye.swap),
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
            lambda.eps1=double(length*n[1]),
            as.double(lambda.eps2),
            lambda.eps2=double(length*n[1]),
            as.double(a.eps),
            as.double(b.eps),
            a.eps=double(length),
            b.eps=double(length),
            as.double(w.in),
            as.double(df.choice),
            as.integer(length(df.choice)),
            as.double(df.in),
            df=double(n[2]*length),
            w=as.double(w.out),
            as.double(rho),
            rho=double(length),
            as.double(lambda.gamma1),
            as.double(lambda.gamma2),
            lambda.gamma1=double(length),
            lambda.gamma2=double(length),
            as.integer(min.iter),
            as.integer(batch),
            as.integer(all.out))


### Create a new object
    new.mcmc<-list(gamma1=t(matrix(obj$gamma1,length,n[1],byrow=TRUE)),
                   gamma2=t(matrix(obj$gamma2,length,n[1],byrow=TRUE)),
                   mu=obj$mu, beta2=obj$beta2, alpha2=obj$alpha2,delta22=obj$delta22,
                   lambda.eps1=t(matrix(obj$lambda.eps1,length,n[1],byrow=TRUE)),
                   lambda.eps2=t(matrix(obj$lambda.eps2,length,n[1],byrow=TRUE)),
                   lambda.gamma1=obj$lambda.gamma1,lambda.gamma2=obj$lambda.gamma2,rho=obj$rho,
                   w=matrix(obj$w,n[1],n[2]),
                   df=t(matrix(obj$df,length,n[2],byrow=TRUE)),
                   eta=t(matrix(obj$eta,length,n[2],byrow=TRUE)),
                   a.eps=obj$a.eps,b.eps=obj$b.eps,shift=shift)
### Give it the right class
    class(new.mcmc)<-"mcmc"
    
    return(new.mcmc)
}

