\name{arrange.row}
\alias{arrange.row}
\title{
  Reorder a dataset by increasing row order.}

\usage{arrange.row(data)}
\description{The functions could be used to reorder a dataset to make
  sure that all the genes are in the same row before fitting any model.
  The \code{arrange.row} function is also used by the \code{weight.plot}
  function to map all the genes to their position on the slide.
}
\arguments{
  \item{data}{
    A dataset containing the row indices in the first column and the
    column indices in the second column.
    The row indices should all be distinct. All indices should start at zero!
  }  
}

\value{ The ordered dataset.}

\references{Robust Estimation of cDNA Microarray Intensities with Replicates
Raphael Gottardo, Adrian E. Raftery, Ka Yee Yeung, and Roger Bumgarner
Department of Statistics, University of Washington, Box 354322, Seattle, WA 98195-4322} 

\seealso{
  \code{weight.plot}
}
\examples{
data(hiv)
### Put the indices in the first two columns and
### reorder the first 4 replicates
new.data<-cbind(hiv[,9:10],hiv[,1:4])
ordered.data<-arrange.row(new.data)
}

\author{Raphael Gottardo}

\keyword{data}
