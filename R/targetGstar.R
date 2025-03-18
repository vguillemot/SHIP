#' Computation of the target Gstar.
#' 
#' The \eqn{p \times p}{p x p} target Gstar is computed from the \eqn{n \times
#' p}{n x p} data matrix. It it a modified version of target G. In particular,
#' it involves two parameters for the correlation (a positive and a negative
#' one) instead of the single parameter \eqn{\bar{r}}{r} in order to account
#' for negatively correlated genes within the same pathway
#' 
#' 
#' @param x A \eqn{n \times p} data matrix.
#' @param genegroups A list of genes obtained using the database KEGG, where
#' each entry itself is a list of pathway names this genes belongs to. If a
#' gene does not belong to any gene functional group, the entry is NA.
#' @return A \eqn{p \times p}{p x p} matrix.
#' @author Monika Jelizarow and Vincent Guillemot
#' @seealso \code{\link{targetCor}}, \code{\link{targetF}},
#' \code{\link{targetG}}, \code{\link{targetGstar}}, \code{\link{targetGpos}}.
#' @references \itemize{ \item J. Schaefer and K. Strimmer, 2005. A shrinkage
#' approach to large-scale covariance matrix estimation and implications for
#' functional genomics.  Statist. Appl. Genet. Mol. Biol. 4:32.  \item M.
#' Jelizarow, V. Guillemot, A. Tenenhaus, K. Strimmer, A.-L. Boulesteix, 2010.
#' Over-optimism in bioinformatics: an illustration. Bioinformatics. Accepted.
#' }
#' @keywords methods multivariate
#' @examples
#' 
#' # A short example on a toy dataset
#' # require(SHIP)
#' data(expl)
#' attach(expl)
#' tar <- targetGstar(x,genegroups)
#' which(tar[upper.tri(tar)]!=0) # not many non zero coefficients !
#' 
#' @importFrom stats sd cov cor
#' @export
targetGstar <- function(x,genegroups) {
  T1   <- target.help(genegroups)
  T2   <- matrix(nrow=length(genegroups),ncol=length(genegroups),data=0)
  corm <- T1*cor(x)
  diag(corm)<-0
  sum.pos <- sum(corm[T1==1 & corm > 0])
  sum.neg <- sum(corm[T1==1 & corm < 0])
  
  cora.pos <- sum.pos/sum(corm>0)
  cora.neg <- sum.neg/sum(corm<0)
  
  for (i in 1:length(genegroups)) {
    for (j in 1:i) {
      if (i!=j & T1[i,j]==1 & corm[i,j] > 0) T2[i,j] <- cora.pos*(sd(x[,i])*sd(x[,j]))
      if (i!=j & T1[i,j]==1 & corm[i,j] < 0) T2[i,j] <- cora.neg*(sd(x[,i])*sd(x[,j]))
      if (i==j)T2[i,j] <- cov(x[,i],x[,j])
      T2[j,i]<-T2[i,j]
    }
  }
  T2
}

