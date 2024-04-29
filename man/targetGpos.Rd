\name{targetGpos}
\alias{targetGpos}
\title{
Computation of the target Gpos.
}
\description{
The \eqn{p \times p}{p x p} target Gpos is computed from the 
\eqn{n \times p}{n x p} data matrix. It it a modified version 
of target G. In particular, it completely ignores negative correlations and 
computes the mean correlation \eqn{\bar{r}}{r} using the positive ones only.
}
\usage{
targetGpos(x, genegroups)
}

\arguments{
  \item{x}{A \eqn{n \times p}{n x p} data matrix.
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
\itemize{
\item J. Schaefer and K. Strimmer, 2005. A shrinkage approach to large-scale 
covariance matrix estimation and implications for functional genomics. 
Statist. Appl. Genet. Mol. Biol. 4:32.
\item M. Jelizarow, V. Guillemot, A. Tenenhaus, K. Strimmer, A.-L. Boulesteix, 2010. 
 Over-optimism in bioinformatics: an illustration. Bioinformatics. Accepted.
}
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
tar <- targetGpos(x,genegroups)
which(tar[upper.tri(tar)]!=0) # not many non zero coefficients !
}

\keyword{methods}
\keyword{multivariate}
