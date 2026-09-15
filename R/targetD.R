#' Diagonal covariance target.
#' 
#' The target keeps the sample variances and sets all covariances to zero. If
#' \eqn{S = (s_{ij})}{S = (s_ij)} is the sample covariance matrix, the target
#' is
#' \deqn{T_{ij} = \begin{cases} s_{ii} & \text{if } i = j \\ 0 & \text{otherwise} \end{cases}}{T_ij = s_ii if i = j, and 0 otherwise.}
#' 
#' @param x A \eqn{n \times p}{n x p} data matrix.
#' @param genegroups The genegroups are not used for this target.
#' @return A \eqn{p \times p}{p x p} diagonal matrix.
#' @author Monika Jelizarow and Vincent Guillemot
#' @family covariance targets
#' @keywords methods multivariate
#' @examples
#' 
#' x <- matrix(rnorm(10*30),10,30)
#' target_matrix <- targetD(x, NULL)
#' 
#' @importFrom stats var
#' @export
targetD <- function(x, genegroups = NULL) {
	validate_data(x)
	diag(apply(x, 2L, var))
}
