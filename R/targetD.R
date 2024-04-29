#' Computation of the diagonal target D ('diagonal, unequal variances').
#' 
#' The \eqn{p \times p}{n x p} diagonal target D is computed from the \eqn{n
#' \times p}{n x p} data matrix. It is defined as follows (\eqn{i,j =
#' 1,...,p}{i,j=1,...,p}): \deqn{t_{ij}=\left\{ }{tij = sii if i=j, 0
#' otherwise}\deqn{\begin {array} {ll} }{tij = sii if i=j, 0
#' otherwise}\deqn{s_{ii}\;&\mbox{if}\;i=j\\ }{tij = sii if i=j, 0
#' otherwise}\deqn{0\;&\mbox{if}\;i\neq j\\ }{tij = sii if i=j, 0
#' otherwise}\deqn{\end {array} }{tij = sii if i=j, 0
#' otherwise}\deqn{\right.}{tij = sii if i=j, 0 otherwise} where
#' \eqn{s_{ij}}{sij} denotes the entry of the unbiased covariance matrix in row
#' \eqn{i}{i}, column \eqn{j}{j}.
#' 
#' 
#' @param x A \eqn{n \times p}{n x p} data matrix.
#' @param genegroups The genegroups are not used for this target.
#' @return A \eqn{p \times p}{p x p} diagonal matrix.
#' @author Monika Jelizarow and Vincent Guillemot
#' @seealso \code{\link{targetCor}}, \code{\link{targetF}},
#' \code{\link{targetG}}, \code{\link{targetGstar}}, \code{\link{targetGpos}}.
#' @references J. Schaefer and K. Strimmer, 2005. A shrinkage approach to
#' large-scale covariance matrix estimation and implications for functional
#' genomics.  Statist. Appl. Genet. Mol. Biol. 4:32.
#' @keywords methods multivariate
#' @examples
#' 
#' x <- matrix(rnorm(10*30),10,30)
#' tar <- targetD(x,NULL)
#' 
targetD <-
function(x,genegroups)  diag(apply(x,2,var))

