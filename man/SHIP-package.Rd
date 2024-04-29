\name{SHIP-package}
\alias{SHIP-package}
\docType{package}
\title{
SHrinkage covariance Incorporating Prior knowledge
}
\description{
The SHIP-package implements the shrinkage estimator of a covariance matrix
given any covariance target, such as described by Schaefer and Strimmer in 2005. 
In addition, it proposes several targets based on biological knowledge extracted 
from the public database KEGG.
}
\details{
To use the shrinkage estimator, one should just have at hand a data set in the form
of a \eqn{n \times p}{n x p} matrix, and a covariance target.

If one wishes to use the proposed targets, the data set should be compatible with KEGG,
i.e. it should be possible to extract for each gene the pathways it belongs to. 
This information, for example, can be found in libraries such as hgu133plus2.db.
}
\author{
Monika Jelizarow and Vincent Guillemot
}
\references{
\itemize{
\item J. Schaefer and K. Strimmer, 2005. A shrinkage approach to large-scale 
covariance matrix estimation and implications for functional genomics. 
Statist. Appl. Genet. Mol. Biol. 4:32.
\item M. Jelizarow, V. Guillemot, A. Tenenhaus, K. Strimmer, A.-L. Boulesteix, 2010. 
 Over-optimism in bioinformatics: an illustration. Bioinformatics. Accepted.
}}
\keyword{ package }
\examples{
# A short example on a toy dataset
# require(SHIP)

data(expl)
attach(expl)

sig1 <- shrink.estim(x,targetD(x))
sig2 <- shrink.estim(x,targetF(x))
sig3 <- shrink.estim(x,targetCor(x,genegroups))
sig4 <- shrink.estim(x,targetG(x,genegroups))

paste(sig1[[2]],collapse=" ")
paste(sig2[[2]],collapse=" ")
paste(sig3[[2]],collapse=" ")
paste(sig4[[2]],collapse=" ")

\dontrun{
# Example on how to get the gene groups lists
require(hgu95av2.db)
# e.g. we have some interesting gene names :
vec <- c("MYC","ID2","PTGER4","ATF4","FGFR1","MET","HLA-DRB6")
# we then want to convert them into Probe Sets
symb <- as.list(hgu95av2SYMBOL)
pbsets <- names(symb[unlist(sapply(vec,function(x,l) which(l==x)[1],symb))])
# Probe Sets which are themselves converted into a gene groups list
genegroups <- as.list(hgu95av2PATH)[pbsets]
}
}