\name{mat.mean}
\alias{mat.mean}
\title{Compute the mean and standard deviation of each row in a data matrix}
\usage{mat.mean(data)
}
\description{
  This function computes the mean and standard deviation of each row in
  a data matrix. The source code is written in C. As a consequence, the
  computation is quite fast.
  
}
\arguments{
  \item{data}{The data matrix.}
}
\value{A matrix, the first columns contain the means, the second the
  standard deviations.
}
%\details{
%}
%\references{
%}
\author{Raphael Gottardo}
\seealso{\code{\link{mean}},\code{\link{sd}}}
\examples{
data(hiv)
sample1<-hiv[,1:4]
ms1<-mat.mean(sample1)
}
\keyword{arith}
