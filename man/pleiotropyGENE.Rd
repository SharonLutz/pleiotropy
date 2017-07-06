\name{pleiotropyGENE}
\alias{pleiotropyGENE}
\title{pleiotropyGENE}
\description{
This function tests for pleiotropy using 2 different gene based approaches
}
\usage{
pleiotropyGENE(X, Y, Ydist, Z = NULL, covariates = FALSE, nPerm = 1000)
}
\arguments{
  \item{X}{
X is a matrix of rare and/or common variants in a region where the number of rows of the matrix equal the number of subjects and the number of columns equal the number of SNPs in the region.
}
  \item{Y}{
Y is a matrix of the traits of interest where the number of rows equal the number of subjects and the number of columns equal the number of traits.
}
  \item{Ydist}{
Ydist is a vector that specifies the distribution of the each trait. For example, for one normally distributed trait and the second binary trait, then Ydist<-c("C","D"). These choices are specified by the SKAT function.
}
  \item{Z}{
Z is a matrix of covariates where the number of rows equal the number of subjects and the number of columns equal the number of covariates.
}
  \item{covariates}{
If covariates=FALSE, then the models will not be adjusted for covariates. If covariates=TRUE, then the model will be adjusted for covariates.
}
  \item{nPerm}{
nPerm is the number of permutations. Default is 1,000.
}
}
\details{
 For gene based tests of pleiotropy, this function computes the p-values obtained from both the cut-off based permutation approach and the Hausdorff based permutation approach.
}
\value{
cutoffPvalue is the p-value from the cut-off based permutation approach

hausdorffPvalue is the p-value from the Hausdorff based permutation approach
}
\author{
Sharon Lutz
}
\section{Warning}{
	SKAT R package is needed to run this function. Make sure the SKAT package has been installed first.
}

\examples{
set.seed(1)
X1<-rbinom(1000,2,0.005)
X2<-rbinom(1000,2,0.005)
X3<-rbinom(1000,2,0.005)
X4<-rbinom(1000,2,0.005)
X5<-rbinom(1000,2,0.005)
X<-cbind(X1,X2,X3,X4,X5)
Z<-matrix(rnorm(1000),nrow=1000,ncol=1)
Y1<-rnorm(1000,0.45*(X1+X2+X3+X4+X5)+0.1*Z,1)
Y2<-rnorm(1000,0.45*(X1+X2+X3+X4+X5)+0.1*Z,1)
Y<-cbind(Y1,Y2)
Ydist<-c("C","C")

pleiotropyGENE(X,Y,Ydist) # NOT adjusting for covariates Z
pleiotropyGENE(X,Y,Ydist,Z,covariates=TRUE) # Adjusting for covariates Z
}

\keyword{ function }
\keyword{ pleiotropy }
