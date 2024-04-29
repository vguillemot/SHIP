\name{targetCor}
\alias{targetCor}
\title{
Computation of the target Cor.
}
\description{
The \eqn{p \times p}{p x p} target Cor is computed from the 
\eqn{n \times p}{n x p} data matrix. It it a modified version of 
target G. In particular, it tests the correlations (with a significance level of 
0.05) and sets the non-significant correlations to zero before the mean correlation 
\eqn{\bar{r}}{r} is computed.
}
\usage{
targetCor(x, genegroups)
}

\arguments{
  \item{x}{A \eqn{n \times p} data matrix.
}
  \item{genegroups}{A list of genes obtained using the database KEGG, where each 
  entry itself is a list of pathway names this genes belongs to. If a gene does 
  not belong to any gene functional group, the entry is NA.
}
}

\value{
A \eqn{p \times p}{p x p} matrix.
}

\references{
J. Schaefer and K. Strimmer, 2005. A shrinkage approach to large-scale 
covariance matrix estimation and implications for functional genomics. 
Statist. Appl. Genet. Mol. Biol. 4:32.
}
\author{
Monika Jelizarow and Vincent Guillemot
}

\seealso{
\code{\link{targetCor}}, 
\code{\link{targetF}}, 
\code{\link{targetG}}, 
\code{\link{targetGstar}}, 
\code{\link{targetGpos}}. 
}

\examples{
# A short example on a toy dataset
# require(SHIP)
data(expl)
attach(expl)
tar <- targetCor(x,genegroups)
which(tar[upper.tri(tar)]!=0) # not many non zero coefficients !
}

\keyword{methods}
\keyword{multivariate}

