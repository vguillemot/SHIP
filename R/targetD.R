#' Computation of the diagonal target D ('diagonal, unequal variances').
#' 
#' The \eqn{p \times p}{n x p} diagonal target D is computed from the \eqn{n
#' \times p}{n x p} data matrix. It is defined as follows (\eqn{i,j =
#' 1,...,p}{i,j=1,...,p}): 
#' \deqn{t_{ij}=\begin{cases}s_{ii} & \text{ if } i=j \\ 0 & \text{ otherwise }\end{cases}}
#' where
#' \eqn{s_{ij}}{sij} denotes the entry of the unbiased covariance matrix in row
#' \eqn{i}{i}, column \eqn{j}{j}.
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
#' @importFrom stats var
#' @export
targetD <- function(x,genegroups)  diag(apply(x,2,var))
