#' Shrinkage estimator of the covariance matrix, given a data set and a
#' covariance target.
#' 
#' The estimator combines the sample covariance matrix \eqn{S}{S} with a
#' target matrix \eqn{T}{T}:
#' \deqn{\widehat{\Sigma} = \lambda T + (1 - \lambda)S}{Sigma_hat = lambda * T + (1 - lambda) * S}
#' The estimated intensity \eqn{\lambda}{lambda} is constrained to
#' \eqn{[0, 1]}{[0, 1]}. Values near zero keep more of the sample covariance;
#' values near one give more weight to the target.
#' 
#' 
#' @param x A \eqn{n \times p}{n x p} matrix (the data set) .
#' @param target A \eqn{p \times p}{p x p} matrix (the covariance target).
#' @return A list with components `shrink.cov`, the estimated covariance
#' matrix, and `lambda`, the estimated shrinkage intensity.
#' @author Monika Jelizarow and Vincent Guillemot
#' @references J. Schaefer and K. Strimmer, 2005. A shrinkage approach to
#' large-scale covariance matrix estimation and implications for functional
#' genomics.  Statist. Appl. Genet. Mol. Biol. 4:32.
#' Jelizarow, M., Guillemot, V., Tenenhaus, A., Strimmer, K. and Boulesteix,
#' A.-L. (2010). Over-optimism in bioinformatics: an illustration.
#' Bioinformatics.
#' @family main functions
#' @keywords methods multivariate
#' @examples
#' 
#' # Simulate dataset
#' x <- matrix(rnorm(20*30),20,30)
#' # Try different targets
#' shrink.estim(x, target = build.target(x, type = "D"))
#' 
#' @importFrom stats cov cor cov2cor
#' @export

shrink.estim <- function(x, target) {
    validate_data(x)
    validate_target(target, ncol(x))

    n_variables <- ncol(x)
    n_observations <- nrow(x)
    covariance_matrix <- cov(x)
    correlation_matrix <- cov2cor(covariance_matrix)
    scaled_data <- scale(x, center = TRUE, scale = TRUE)
    scaled_crossprod <- crossprod(scaled_data)
    squared_scaled_crossprod <- scaled_crossprod^2

    variance_term <- (n_observations / ((n_observations - 1)^3)) *
        (crossprod(scaled_data^2) - 1 / n_observations * squared_scaled_crossprod)
    diag(variance_term) <- 0

    moment_matrix <- matrix(
        rep(apply(scaled_data^2, 2, mean), n_variables),
        n_variables,
        n_variables
    )
    fourth_moment_term <- (n_observations / (2 * (n_observations - 1)^3)) *
        (crossprod(scaled_data^3, scaled_data) +
            crossprod(scaled_data, scaled_data^3) -
            (moment_matrix + t(moment_matrix)) * scaled_crossprod)
    diag(fourth_moment_term) <- 0
    fourth_moment_term[target == 0] <- 0

    target_correlation <- cov2cor(target)
    squared_difference <- (correlation_matrix - target_correlation)^2
    denominator <- sum(squared_difference)
    lambda <- if (denominator <= .Machine$double.eps) {
        0
    } else {
        (sum(variance_term) - sum(target_correlation * fourth_moment_term)) /
            denominator
    }
    lambda <- max(min(lambda, 1), 0)
    shrinkage_covariance <- lambda * target + (1 - lambda) * covariance_matrix

    list(shrink.cov = shrinkage_covariance, lambda = lambda)
}

