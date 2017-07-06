\name{hausdorff}
\alias{hausdorff}
\title{hausdorff}
\description{
This function computes the Hausdorff metric
}
\usage{
hausdorff(Yob, Yperm)
}
\arguments{
  \item{Yob}{
  	Let Yob be a vector of numbers (i.e. p-values for the observed traits with the SNP of interest)
}
  \item{Yperm}{
    Let Yperm be a vector of numbers (i.e. p-values for the permuted traits with the SNP of interest)
}
}
\details{
 This function compares the observed set of p-values to the permuted using the hausdorff function, max(sup_x inf_y d(x,y), sup_y inf_x d(x,y))
}
\value{
hausdorff metric for the 2 sets
}

\author{
Sharon Lutz
}

\examples{
hausdorff(c(0.05,0.01),c(0.001,0.9))
}

\keyword{ function }
\keyword{ hausdorff }
