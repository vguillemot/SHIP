#' Constant-correlation covariance target.
#' 
#' The target keeps the sample variances and replaces every off-diagonal
#' correlation by their average. If \eqn{S = (s_{ij})}{S = (s_ij)} is the
#' sample covariance matrix and \eqn{\bar{r}}{r_bar} is the mean sample
#' correlation, then
#' \deqn{
#'   T_{ij} = \begin{cases}
#'    s_{ii} & \text{if } i = j \\
#'    \bar{r}\sqrt{s_{ii}s_{jj}} & \text{otherwise}
#'   \end{cases}
#' }{
#'   T_ij = s_ii if i = j, and r_bar * sqrt(s_ii * s_jj) otherwise.
#' }
#' 
#' @param x A \eqn{n \times p}{n x p} data matrix.
#' @param genegroups The genegroups are not used for this target.
#' @return A \eqn{p \times p}{p x p} matrix.
#' @author Monika Jelizarow and Vincent Guillemot
#' @family covariance targets
#' @keywords methods multivariate
#' @examples
#' 
#' # A short example on a toy dataset
#' # require(SHIP)
#' data(expl)
#' target_matrix <- targetF(expl$x, NULL)
#' which(target_matrix[upper.tri(target_matrix)] != 0) # many non zero coefficients !
#' 
#' @importFrom stats cor cov
#' @export
targetF <- function(x, genegroups = NULL) {
  validate_data(x)
  covariance_matrix <- cov(x)
  correlation_matrix <- cor(x)
  correlations <- correlation_matrix[upper.tri(correlation_matrix)]
  mean_correlation <- if (length(correlations) == 0L) 0 else mean(correlations)
  target_matrix <- mean_correlation * tcrossprod(sqrt(diag(covariance_matrix)))
  diag(target_matrix) <- diag(covariance_matrix)
  target_matrix
}
