#' Computation of target G ('knowledge-based constant correlation model').
#' 
#' The \eqn{p \times p}{p x p} target G is computed from the \eqn{n \times p}{n
#' x p} data matrix. It is defined as follows (\eqn{i,j = 1,...,p}{i,j =
#' 1,...,p}): \deqn{t_{ij}=\left\{ }{sii if i=j ; r\sqrt{s_{ii}s_{jj}} if i~j ;
#' 0 otherwise }\deqn{\begin {array} {ll} }{sii if i=j ; r\sqrt{s_{ii}s_{jj}}
#' if i~j ; 0 otherwise }\deqn{s_{ii}\;&\mbox{if}\;i=j\\ }{sii if i=j ;
#' r\sqrt{s_{ii}s_{jj}} if i~j ; 0 otherwise
#' }\deqn{\bar{r}\sqrt{s_{ii}s_{jj}}\;&\mbox{if}\;i\neq j, i\sim j\\ }{sii if
#' i=j ; r\sqrt{s_{ii}s_{jj}} if i~j ; 0 otherwise }\deqn{0\;&\mbox{otherwise}
#' }{sii if i=j ; r\sqrt{s_{ii}s_{jj}} if i~j ; 0 otherwise }\deqn{\end{array}
#' }{sii if i=j ; r\sqrt{s_{ii}s_{jj}} if i~j ; 0 otherwise }\deqn{\right.}{sii
#' if i=j ; r\sqrt{s_{ii}s_{jj}} if i~j ; 0 otherwise } where \eqn{\bar{r}}{r}
#' is the average of sample correlations and \eqn{s_{ij}}{sij} denotes the
#' entry of the unbiased covariance matrix in row \eqn{i}{i}, column
#' \eqn{j}{j}. The notation \eqn{i\sim j}{i ~ j} means that genes \eqn{i}{i}
#' and \eqn{j}{j} are connected, i.e. genes \eqn{i}{i} and \eqn{j}{j} are in
#' the same gene functional group.
#' 
#' 
#' @param x A \eqn{n \times p}{n x p} data matrix.
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
#' tar <- targetG(x,genegroups)
#' which(tar[upper.tri(tar)]!=0) # not many non zero coefficients !
#' 
targetG <-
function(x,genegroups) {
  T1   <- target.help(genegroups)
  T2   <- matrix(nrow=length(genegroups),ncol=length(genegroups),data=0)
  corm <- T1*cor(x)
  diag(corm)<-0
  cora <-sum(colSums(corm))/sum(corm!=0)
  for (i in 1:length(genegroups)) {
    for (j in 1:i) {
      if (i!=j & T1[i,j]==1) T2[i,j] <- cora*(sd(x[,i])*sd(x[,j]))
      if (i==j)T2[i,j] <- cov(x[,i],x[,j])
      T2[j,i]<-T2[i,j]
      }
  }
  T2
}

