ls.effect<-function(sample1,sample2,dye.swap=FALSE,nb.col1=NULL)
  {
    n<-dim(sample1)
    n1<-n[1]
    n2<-n[2]
### Gene effect sample 1
    gamma1<-rep(0,n1)
### Gene effect sample 2
    gamma2<-rep(0,n1)
    
### Array effect
    eta<-rep(0,n2)

### Residuals
    R1<-sample1
    R2<-sample2
### Main effects
    M1<-sample1
    M2<-sample2

    ### Some constants
    m11<-mean(as.double(sample1[,1:nb.col1]),na.rm=TRUE)
    m12<-mean(as.double(sample1[,(nb.col1+1):n2]),na.rm=TRUE)
    m21<-mean(as.double(sample2[,1:nb.col1]),na.rm=TRUE)
    m22<-mean(as.double(sample2[,(nb.col1+1):n2]),na.rm=TRUE)

    if(dye.swap)
      {
### Offset
        mu<-m11
### Compute the dye effect
        beta2<-m12-m11
        
### Compute the sample effect
        alpha2<-m21-m11
        
### delta22 
        delta22<-m22+m11-m12-m21
        
### Compute the array effect
        
        for(i in 1:nb.col1)
          {
            eta[i]<-1/2*(mean(sample1[,i],na.rm=TRUE)+mean(sample2[,i],na.rm=TRUE)-m11-m21)
          }
        for(i in (nb.col1+1):n2)
          {
            eta[i]<-1/2*(mean(sample1[,i],na.rm=TRUE)+mean(sample2[,i],na.rm=TRUE)-m12-m22)
          }
        
        
### Compute the gene effect
        ms1<-mat.mean(sample1)[,1]
        ms2<-mat.mean(sample2)[,1]
        for(i in 1:n1)
          {
            gamma1[i]<-ms1[i]-(m11+m12)/2.
            gamma2[i]<-ms2[i]-(m21+m22)/2.
            for(j in 1:nb.col1)
              {
                M1[i,j]<-mu+eta[j]+gamma1[i]
                M2[i,j]<-mu+eta[j]+alpha2+gamma2[i]
                
                R1[i,j]<-sample1[i,j]-M1[i,j]
                R2[i,j]<-sample2[i,j]-M2[i,j]
              }
            for(j in (nb.col1+1):n2)
              {
                M1[i,j]<-mu+eta[j]+beta2+gamma1[i]
                M2[i,j]<-mu+eta[j]+alpha2+beta2+delta22+gamma2[i]
                R1[i,j]<-sample1[i,j]-M1[i,j]
                R2[i,j]<-sample2[i,j]-M2[i,j]
              } 
          }     
      }
    else
      {
        beta2<-NULL
        delta22<-NULL
### Offset
        mu<-1/2*(m11+m12) 
        
### Compute the sample effect
        alpha2<-1/2*(m21+m22-m11-m12)
        
### delta22 
        delta22<-1/2*(m22+m11-m12-m21)
        
### Compute the array effect

        for(i in 1:(n2))
          {
            eta[i]<-1/2*(mean(sample1[,i],na.rm=TRUE)+mean(sample2[,i],na.rm=TRUE)-(m11+m12)/2-(m21+m22)/2)
          }
        
        
        
### Compute the gene effect
        ms1<-mat.mean(sample1)[,1]
        ms2<-mat.mean(sample2)[,1]
        for(i in 1:n1)
          {
            gamma1[i]<-ms1[i]-(m11+m12)/2.
            gamma2[i]<-ms2[i]-(m21+m22)/2.
            for(j in 1:nb.col1)
              {
                M1[i,j]<-mu+eta[j]+gamma1[i]
                M2[i,j]<-mu+eta[j]+alpha2+gamma2[i]
                R1[i,j]<-sample1[i,j]-M1[i,j]
                R2[i,j]<-sample2[i,j]-M2[i,j]
              }
            for(j in (nb.col1+1):n2)
              {
                M1[i,j]<-mu+eta[j]+gamma1[i]
                M2[i,j]<-mu+eta[j]+alpha2+gamma2[i]
                R1[i,j]<-sample1[i,j]-M1[i,j]
                R2[i,j]<-sample2[i,j]-M2[i,j]
              }
          }
      }
    
    list(mu=mu,eta=eta,alpha2=alpha2,beta2=beta2,delta22=delta22,gamma1=gamma1,gamma2=gamma2,M1=M1,M2=M2,R1=R1,R2=R2)
  }
