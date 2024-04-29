\name{targetF}
\alias{targetF}
\title{
Computation of target F ('constant correlation model').
}
\description{
The \eqn{p \times p}{p x p} target F is computed from the 
\eqn{n \times p}{n x p} data matrix.
It is defined as follows (\eqn{i,j = 1,...,p}{i,j = 1,...,p}):
\deqn{t_{ij}=\left\{
\begin {array} {ll}
s_{ii}\;&\mbox{if}\;i=j\\
\bar{r}\sqrt{s_{ii}s_{jj}}\;&\mbox{if}\;i\neq j\\
\end {array}
\right.}{tij = sii if i=j ; r\sqrt{s_{ii}s_{jj}} otherwise} where \eqn{\bar{r}}{r} 
is the average of sample correlations and \eqn{s_{ij}}{sij} denotes the entry of 
 the unbiased covariance matrix in row \eqn{i}{i}, column \eqn{j}{j}.
}
\usage{
targetF(x, genegroups)
}
\arguments{
  \item{x}{A \eqn{n \times p}{n x p} data matrix.
}
  \item{genegroups}{The genegroups are not used for this target.
}
}
\value{
A \eqn{p \times p}{p x p} matrix.
}

\references{
J. Schaefer and K. Strimmer, 2005. A shrinkage approach to large-scale 
covariance matrix estimation and implications for functional genomics. 
Statist. Appl. Genet. Mol. Biol. 4:32. }

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
tar <- targetF(x,NULL)
which(tar[upper.tri(tar)]!=0) # many non zero coefficients !
}

\keyword{methods}
\keyword{multivariate}
