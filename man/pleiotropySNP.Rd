\name{pleiotropySNP}
\alias{pleiotropySNP}
\title{pleiotropySNP}
\description{
This function tests for pleiotropy using 2 different SNP based approaches
}
\usage{
pleiotropySNP(X, Y, Ydist, Z = NULL, covariates = FALSE, nPerm = 5000)
}
\arguments{
  \item{X}{
  	X is the vector of the SNP of interest
}
  \item{Y}{
  	Y is a matrix of the traits of interest where the number of rows equal the number of subjects and the number of columns equal the number of traits.
}
  \item{Ydist}{
  	Ydist is a vector that specifies the distribution of the each trait. For example, for one normally distributed trait and the second binary trait, then Ydist<-c("gaussian","binomial"). Other choices for Ydist can be found by looking at the choice of family for the glm function (i.e. ?glm).
}
  \item{Z}{
  	Z is a matrix of covariates where the number of rows equal the number of subjects and the number of columns equal the number of covariates.
}
  \item{covariates}{
  	If covariates=FALSE, then the models will not be adjusted for covariates. If covariates=TRUE, then the model will be adjusted for covariates.
}
  \item{nPerm}{
  	nPerm is the number of permutations. Default is 5,000.
}
}
\details{
For SNP based tests of pleiotropy, this function computes the p-values obtained from both the cut-off based permutation approach and the Hausdorff based permutation approach.
}
\value{
cutoffPvalue is the p-value from the cut-off based permutation approach

hausdorffPvalue is the p-value from the Hausdorff based permutation approach
}
\author{
Sharon Lutz
}
\examples{
set.seed(1)
X<-rbinom(1000,2,0.2) #generate the SNP X=0,1,2
Z<-matrix(rnorm(1000),nrow=1000,ncol=1)
Y1<-rnorm(1000,0.2*X+0.1*Z,1)
Y2<-rnorm(1000,0.2*X+0.1*Z,1)
Y<-cbind(Y1,Y2)
Ydist<-c("gaussian","gaussian")

pleiotropySNP(X,Y,Ydist) # NOT adjusting for covariates Z
pleiotropySNP(X,Y,Ydist,Z,covariates=TRUE) # Adjusting for covariates Z
}
\keyword{ function }
\keyword{ pleiotropy }
