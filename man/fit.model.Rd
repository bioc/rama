\name{fit.model}
\alias{fit.model}
\title{Robust estimation of microarray intensities with replicates}
\usage{fit.model(sample1,sample2,B=1000,min.iter=0,batch=10,shift=NULL,mcmc.obj=NULL,dye.swap=FALSE,nb.col1=NULL,all.out=FALSE,ci=0.95,
verbose=FALSE)
     }
\description{
  Estimate the log transformed intensities of each sample of a replicated
  microarray experiment. The estimation is done via Hiearchical Bayesian
  Modeling.
}
\arguments{
  \item{sample1}{
    The matrix of intensity from the sample 1. Each row corresponds to a different gene.}
  \item{sample2}{
    The matrix of intensity from the sample 2. Each row corresponds to a different gene.}
  \item{B}{
    The number of iteration used the MCMC algorithm.}
  \item{min.iter}{
    The length of the burn-in period in the MCMC algorithm.\code{min.iter} should be less than
    B.
  }
  \item{batch}{The thinning value to be used in the MCMC. Only every
    \code{batch}-th iteration will be stored.}
  \item{mcmc.obj}{An object of type mcmc, as returned by \code{fit.model}. \code{mcmc.obj} is used to initialized the MCMC.
    If no \code{mcmc.obj}, the MCMC is initialized to the least squares estimates.}
  \item{shift}{The shift to be used in the log transformation. If
    \code{shift=NULL} is specified (default), it is estimated using \code{est.shift}}
    \item{dye.swap}{A logical value indicating if the experiment was a dye
    swap experiment.}
  \item{nb.col1}{An integer value corresponding to the number of arrays
    (columns) in the first group of the dye swap experiment. In other words, the number of
    replicates before the dyes have been swaped. }
  \item{all.out}{A logical value indicating if all the parameters
    should be outputted. If \code{all.out} is FALSE, only the posterior
    mean is outputted. This could be used to save memory. }
  \item{ci}{A number between 0 and 1 corresponding to the level used
    when computing log ratio credible intervals. If \code{all.out} is
    FALSE, this option is ignored. }
  \item{verbose}{A logical value indicating if
    the current MCMC iteration number should be printed out.}
}
\value{
  An object of type \code{mcmc} containing the sampled values from the
  posterior distribution.
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
  \item{q.low}{A vector containing the lower bounds for the log ratio
    credible intervals, i.e. the credible intervals for gamma1-gamma2.}
  \item{q.up}{A vector containing the upper bounds for the log ratio
    credible intervals, i.e. the credible intervals for gamma1-gamma2.}
  \item{lambda.gamma1}{
    A vector containing the sampled values for the precision of the gene
    effect prior in sample 1.}
  \item{lambda.gamma2}{
    A vector containing the sampled values for the precision of the gene
    effect prior in sample 2.}
  \item{rho}{A vector containing the sampled values from between sample
    correlation coefficient \code{rho}}
  \item{lambda_eps1}{A matrix, each row contains the sampled values from the
    corresponding gene precision in sample 1.}
  \item{lambda_eps2}{A matrix, each row contains the sampled values from the
    corresponding gene precision in sample 2.}
  \item{a.eps}{A vector containing the sampled values for the mean of
    the prior of the genes precision.}
  \item{b.eps}{A vector containing the sampled values for the variance of
    the prior of the genes precision.}
  \item{w}{A matrix, each element (i,j) correspond to the posterior mean
    of the sampled weights of replicate j in gene i.To save memory, we
    only store the posterior means of the weigths.}
  \item{shift}{The value of the shift.}
  
}

\details{
  The function fits a hierarchical Bayesian model for robust 
  estimation of cDNA microarray intensities.
  Our model addresses classical issues such as design effects, normalization and transformation.
  Outliers are modeled explicitly using a t-distribution.
  Parameter estimation is carried out using Markov Chain Monte Carlo.
}
\references{Robust Estimation of cDNA Microarray Intensities with Replicates
  Raphael Gottardo, Adrian E. Raftery, Ka Yee Yeung, and Roger Bumgarner
  Department of Statistics, University of Washington, Box 354322, Seattle, WA 98195-4322} 
}
\seealso{
  \code{est.shift}
}
\examples{
data(hiv)
mcmc.hiv<-fit.model(hiv[1:10,c(1:4)],hiv[1:10,c(5:8)],B=2000,min.iter=000,batch=1,shift=30,mcmc.obj=NULL,dye.swap=TRUE,nb.col1=2)
}

\author{Raphael Gottardo}

\keyword{models}
\keyword{robust}
