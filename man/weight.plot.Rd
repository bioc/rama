\name{weight.plot}
\alias{weight.plot}
\title{Plot the weights of a given array using the spatial location of
  the genes on the slide}
\usage{weight.plot(mcmc.obj,coordinate,array=1)}
\description{Plot the weights of a given array using the spatial location of
  the genes on the slide.This function is a useful diagnostic tool.
}
\arguments{
  \item{mcmc.obj}{
    An object of class \code{mcmc} as returned by \code{fit.model}}
  \item{coordinate}{
    The coordinate of each gene on the corresponding array. The
    coordinates should be a two column integer valued matrix containing
    the row indices (column 1) and the column indices (column 2).
    The row indices should all be distinct. All indices should start at zero!
  }
  \item{array}{
    An integer corresponding to the array number to be plotted.}
}

\value{ The image plot of the weights. A small weight (bright color)
  correspond to an outlier.
}

\references{Robust Estimation of cDNA Microarray Intensities with Replicates
Raphael Gottardo, Adrian E. Raftery, Ka Yee Yeung, and Roger Bumgarner
Department of Statistics, University of Washington, Box 354322, Seattle, WA 98195-4322} 

\seealso{
  \code{arrange.row}
}
\examples{
data(hiv)
### Initialize the proposals
mcmc.hiv<-fit.model(hiv[1:640,c(1:4)],hiv[1:640,c(5:8)],B=1000,min.iter=500,batch=1,shift=30,mcmc.obj=NULL,dye.swap=TRUE,nb.col1=2)
weight.plot(mcmc.hiv,hiv[1:640,9:10],array=3)
}

\author{Raphael Gottardo}

\keyword{hplot}
