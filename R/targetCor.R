#' Computation of the target Cor.
#' 
#' The \eqn{p \times p}{p x p} target Cor is computed from the \eqn{n \times
#' p}{n x p} data matrix. It it a modified version of target G. In particular,
#' it tests the correlations (with a significance level of 0.05) and sets the
#' non-significant correlations to zero before the mean correlation
#' \eqn{\bar{r}}{r} is computed.
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
#' @references J. Schaefer and K. Strimmer, 2005. A shrinkage approach to
#' large-scale covariance matrix estimation and implications for functional
#' genomics.  Statist. Appl. Genet. Mol. Biol. 4:32.
#' @keywords methods multivariate
#' @examples
#' 
#' # A short example on a toy dataset
#' # require(SHIP)
#' data(expl)
#' attach(expl)
#' tar <- targetCor(x,genegroups)
#' which(tar[upper.tri(tar)]!=0) # not many non zero coefficients !
#' @importFrom stats sd cov cor cor.test
#' @export
targetCor <- function(x,genegroups) {
  
  T1 <- target.help(genegroups)
  T2 <- matrix(nrow=length(genegroups),ncol=length(genegroups),data=0)     
  
  for (i in 2:length(genegroups)) {
    for (j in 1:(i-1)) if (T1[i,j]==1) T1[i,j] <- ifelse(stats::cor.test(x[,i],x[,j])$p.value < 0.05,1,0)
  }
  
  corm <- T1 * cor(x)     
  diag(corm) <- 0  
  
  cora <- ifelse(sum(corm!=0)==0,0,sum(colSums(corm))/sum(corm!=0))
  
  for (i in 1:length(genegroups)) {
    for (j in 1:i) {
      if (i!=j & T1[i,j]==1) T2[i,j] <- cora*(sd(x[,i])*sd(x[,j]))
      if (i==j)T2[i,j] <- cov(x[,i],x[,j])
      T2[j,i]<-T2[i,j]
    }
  }
  T2
}
