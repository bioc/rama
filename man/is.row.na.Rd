\name{is.row.na}
\alias{is.row.na}
\title{Test if a matrix contains missing values}
\usage{is.row.na(data)
     }
\description{
The function returns a vector of logical variables, one for each row of
the matrix. The variable is \code{TRUE} if the row does not contain any
missing values and \code{FAlSE} otherwise.
}
\arguments{
  \item{data}{
    The data matrix.
  }
}
\value{
  The vector of logical variable    
}

\seealso{
  is.na
}
\examples{

### Generate a matrix
M<-matrix(rnorm(100),10,10)
M[1,1]<-NA
M[1,2]<- -Inf
M[3,10]<-NA

### Indices of the rows without missing values 
ind<-is.row.na(M)

### Submatrix of M with finite values
M.finite<-M[ind,]
}

\author{Raphael Gottardo}

\keyword{}
