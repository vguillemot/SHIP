\name{targetG}
\alias{targetG}
\title{
Computation of target G ('knowledge-based constant correlation model').
}
\description{
The \eqn{p \times p}{p x p} target G is computed from the 
\eqn{n \times p}{n x p} data matrix.
It is defined as follows (\eqn{i,j = 1,...,p}{i,j = 1,...,p}):
\deqn{t_{ij}=\left\{
\begin {array} {ll}
s_{ii}\;&\mbox{if}\;i=j\\
\bar{r}\sqrt{s_{ii}s_{jj}}\;&\mbox{if}\;i\neq j, i\sim j\\
0\;&\mbox{otherwise}
\end{array}
\right.}{sii if i=j ; r\sqrt{s_{ii}s_{jj}} if i~j ; 0 otherwise } where \eqn{\bar{r}}{r} 
is the average of sample correlations and \eqn{s_{ij}}{sij} denotes the entry of 
 the unbiased covariance matrix in row \eqn{i}{i}, column \eqn{j}{j}. The 
 notation \eqn{i\sim j}{i ~ j} means that genes
 \eqn{i}{i} and \eqn{j}{j} are connected, i.e. genes \eqn{i}{i} and \eqn{j}{j} 
 are in the same gene functional group.
}

\usage{
targetG(x, genegroups)
}

\arguments{
  \item{x}{A \eqn{n \times p}{n x p} data matrix.
}
  \item{genegroups}{A list of genes obtained using the database KEGG, where each entry itself is a list of pathway names
   this genes belongs to. If a gene does not belong to any gene functional group, the entry is NA.
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
tar <- targetG(x,genegroups)
which(tar[upper.tri(tar)]!=0) # not many non zero coefficients !
}

\keyword{methods}
\keyword{multivariate}