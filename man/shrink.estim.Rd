\name{shrink.estim}
\alias{shrink.estim}
\title{
Shrinkage estimator of the covariance matrix, given a data set and a covariance target.
}
\description{
The shrinkage estimator is computed independently of the target's nature.
}
\usage{
shrink.estim(x, tar)
}
\arguments{
  \item{x}{ A \eqn{n \times p}{n x p} matrix (the data set) .
}
  \item{tar}{A \eqn{p \times p}{p x p} matrix (the covariance target).}
}
\value{
A \eqn{p \times p}{p x p} shrinkage covariance matrix and the estimated \eqn{\lambda}{lambda}.
}
\references{
J. Schaefer and K. Strimmer, 2005. A shrinkage approach to large-scale 
covariance matrix estimation and implications for functional genomics. 
Statist. Appl. Genet. Mol. Biol. 4:32. }
\author{
Monika Jelizarow and Vincent Guillemot
}
\examples{
# Simulate dataset
x <- matrix(rnorm(20*30),20,30)
# Try different targets
shrink.estim(x,tar=build.target(x,type="D"))
}
\keyword{methods}
\keyword{multivariate}
