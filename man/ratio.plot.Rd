\name{ratio.plot}
\alias{ratio.plot}
\title{Plot the estimated log ratios against the overall intensities}
\usage{ratio.plot(mcmc.obj,col=1,pch=1)}
\description{
  Plot the estimated log2(gamma1/gamma2) against log2(gamma1*gamma2)/2.
}
\arguments{
  \item{mcmc.obj}{
    An object of class \code{mcmc} as returned by \code{fit.model}}
  \item{col}{
    The color to be used for the symbols}
  \item{pch}{
    The type of symbols to be used.}
}

\value{ The graph!
}

}
\references{Robust Estimation of cDNA Microarray Intensities with Replicates
Raphael Gottardo, Adrian E. Raftery, Ka Yee Yeung, and Roger Bumgarner
Department of Statistics, University of Washington, Box 354322, Seattle, WA 98195-4322} 
}
\seealso{
  \code{fit.model}
}
\examples{
data(hiv)
### Initialize the proposals
mcmc.hiv<-fit.model(hiv[1:10,c(1:4)],hiv[1:10,c(5:8)],B=2000,min.iter=000,batch=1,shift=30,mcmc.obj=NULL,dye.swap=TRUE,nb.col1=2)
ratio.plot(mcmc.hiv,col=1,pch=1)
}

\author{Raphael Gottardo}

\keyword{hplot}
