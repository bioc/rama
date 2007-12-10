\name{est.shift}
\alias{est.shift}
\title{Estimate the shift used in the log transformation}
\usage{est.shift(sample1,sample2,B=1000,min.iter=0,batch=10,mcmc.obj=NULL,dye.swap=FALSE,nb.col1=NULL,all.out=TRUE,verbose=FALSE)}
\description{
  Estimate the shift in the log transformation when fitting the
  Hierarchical model as in bayes.rob.}
\arguments{
  \item{sample1}{
    The matrix of intensity from the sample 1. Each row corresponds to a different gene.}
  \item{sample2}{
    The matrix of intensity from the sample 2. Each row corresponds to a different gene.}
  \item{B}{
    The number of iteration used the MCMC algorithm.}
  \item{min.iter}{
    The length of the burn-in period in the MCMC algorithm.\code{min.iter} should be less than
    B.}
  \item{batch}{The thinning value to be used in the MCMC. Only every \code{batch}-th iteration will be stored.}
  \item{mcmc.obj}{An object of type \code{mcmc.shift}, as returned by
    \code{est.shift}.
    If no \code{mcmc.obj}, the MCMC is initialized to the least squares estimates.}
  \item{dye.swap}{A logical value indicating if the experiment was a dye
    swap experiment.}
  \item{nb.col1}{An integer value correspinding to the number of arrays
    (columns) in the first group of the dye swap experiment. In other words, the number of
    replicates before the dyes have been swaped. }
  \item{all.out}{A logical value indicating if all the parameters
    should be outputted. If \code{all.out} is FALSE, only the posterior
    mean is outputted. This could be used to save memory. }
  \item{verbose}{A logical value indicating if
    the current MCMC iteration number should be printed out.}
}

\value{
    An object of type \code{mcmc.est} containing the sampled values from the posterior distribution.
  \item{mu}{A vector containing the sampled values from \code{mu}, the baseline intensity.}
  \item{alpha2}{A vector containing the sampled values from \code{alpha2}, the sample effect.}
  \item{beta2}{A vector containing the sampled values from \code{beta2},
    the dye effect.}
  \item{delta22}{A vector containing the sampled values from
    \code{delta_22}, the dye*sample interaction.}
  \item{eta}{A matrix, each row contains the sampled values from the
    corresponding array effect.}
  \item{gamma1}{A matrix, each row contains the sampled values from the
    corresponding gene effect in sample 1.}
  \item{gamma2}{A matrix, each row contains the sampled values from the
    corresponding gene effect in sample 1.}
  \item{lambda.gamma1}{
    A vector containing the sampled values for the precision of the gene
    effect prior in sample 1.}
  \item{lambda.gamma2}{
    A vector containing the sampled values for the precision of the gene
    effect prior in sample 2.}
  \item{rho}{A vector containing the sampled values from between sample
    correlation coefficient \code{rho}}
  \item{lambda_eps1}{A vector containing the sampled values from the
    gene precision in sample 1.}
  \item{lambda_eps2}{A vector containing the sampled values from the
    gene precision in sample 2.}
  \item{shift}{A vector containing the sampled values from the shift.}
}

\details{
  The estimation is done by fitting the same model (as in fit.model) with constant
  variance, Gaussian errors and a prior for the shift.
  The main purpose of this function is to estimate the shift in the log
  transformation.
  Parameter estimation is carried out using Markov Chain Monte
  Carlo. The shift is estimated with the posterior mean.}
\references{Robust Estimation of cDNA Microarray Intensities with Replicates
Raphael Gottardo, Adrian E. Raftery, Ka Yee Yeung, and Roger Bumgarner
Department of Statistics, University of Washington, Box 354322, Seattle, WA 98195-4322} 

\seealso{
\code{fit.model}
}
\examples{
data(hiv)
### Initialize the proposals
mcmc.hiv<-est.shift(hiv[1:10,c(1:4)],hiv[1:10,c(5:8)],B=2000,min.iter=000,batch=1,mcmc.obj=NULL,dye.swap=TRUE,nb.col1=2)
}

\author{Raphael Gottardo}

\keyword{models}

