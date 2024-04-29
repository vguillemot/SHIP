#' Computation of target F ('constant correlation model').
#' 
#' The \eqn{p \times p}{p x p} target F is computed from the \eqn{n \times p}{n
#' x p} data matrix. It is defined as follows (\eqn{i,j = 1,...,p}{i,j =
#' 1,...,p}): \deqn{t_{ij}=\left\{ }{tij = sii if i=j ; r\sqrt{s_{ii}s_{jj}}
#' otherwise}\deqn{\begin {array} {ll} }{tij = sii if i=j ;
#' r\sqrt{s_{ii}s_{jj}} otherwise}\deqn{s_{ii}\;&\mbox{if}\;i=j\\ }{tij = sii
#' if i=j ; r\sqrt{s_{ii}s_{jj}}
#' otherwise}\deqn{\bar{r}\sqrt{s_{ii}s_{jj}}\;&\mbox{if}\;i\neq j\\ }{tij =
#' sii if i=j ; r\sqrt{s_{ii}s_{jj}} otherwise}\deqn{\end {array} }{tij = sii
#' if i=j ; r\sqrt{s_{ii}s_{jj}} otherwise}\deqn{\right.}{tij = sii if i=j ;
#' r\sqrt{s_{ii}s_{jj}} otherwise} where \eqn{\bar{r}}{r} is the average of
#' sample correlations and \eqn{s_{ij}}{sij} denotes the entry of the unbiased
#' covariance matrix in row \eqn{i}{i}, column \eqn{j}{j}.
#' 
#' 
#' @param x A \eqn{n \times p}{n x p} data matrix.
#' @param genegroups The genegroups are not used for this target.
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
#' tar <- targetF(x,NULL)
#' which(tar[upper.tri(tar)]!=0) # many non zero coefficients !
#' 
targetF <-
function(x,genegroups) {
  covm <- cov(x)
  corm <- cor(x)
  var.only <- covm/corm
  diag(corm) <- 0
  cora  <- sum(corm)/sum(corm!=0)
  T <- cora*var.only
  diag(T)<- diag(covm)
  T
}

