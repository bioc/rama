\name{ls.effect}
\alias{ls.effect}
\title{Compute the least squares estimates of the all the effects of the
  general model.}
\usage{
ls.effect(sample1,sample2,dye.swap=FALSE,nb.col1=NULL)
     }
\description{
  Compute the least squares estimates of the all the effects of the
  general model.
}
\arguments{
  \item{sample1}{ 
    The matrix of intensity from the sample 1. Each row corresponds to a different gene.}
  \item{sample2}{
    The matrix of intensity from the sample 2. Each row corresponds to a different gene.}
  \item{dye.swap}{A logical value indicating if the experiment was a dye
    swap experiment.}
  \item{nb.col1}{An integer value correspinding to the number of arrays
    (columns) in the first group of the dye swap experiment. In other words, the number of
    replicates before the dyes have been swaped. } 
}
\value{
  \item{mu}{The baseline intensity}
  \item{alpha2}{The sample effect}
  \item{beta2}{The dye effect}
  \item{delta22}{The dye*sample interaction}
  \item{eta}{The array effects}
  \item{gamma1}{The genes effects in sample 1}
  \item{gamma2}{The genes effect in sample 2}
  \item{M1}{The main effects in sample 1}
  \item{M2}{The main effects in sample 2}
  \item{R1}{The residuals from the sample 1}
  \item{R2}{The residuals from the sample 2}
  
}

\details{
 }
\references{Robust Estimation of cDNA Microarray Intensities with Replicates
  Raphael Gottardo, Adrian E. Raftery, Ka Yee Yeung, and Roger Bumgarner
  Department of Statistics, University of Washington, Box 354322, Seattle, WA 98195-4322} 
}
\seealso{
  \code{fit.model}
}
\examples{
### Compute the least squares effects on the log scale
data(hiv)
ls.fx<-ls.effect(log2(hiv[,c(1:4)]),log2(hiv[,c(5:8)]),dye.swap=TRUE,nb.col1=2)
}

\author{Raphael Gottardo}

\keyword{models}
