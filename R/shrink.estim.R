#' Shrinkage estimator of the covariance matrix, given a data set and a
#' covariance target.
#' 
#' The shrinkage estimator is computed independently of the target's nature.
#' 
#' 
#' @param x A \eqn{n \times p}{n x p} matrix (the data set) .
#' @param tar A \eqn{p \times p}{p x p} matrix (the covariance target).
#' @return A \eqn{p \times p}{p x p} shrinkage covariance matrix and the
#' estimated \eqn{\lambda}{lambda}.
#' @author Monika Jelizarow and Vincent Guillemot
#' @references J. Schaefer and K. Strimmer, 2005. A shrinkage approach to
#' large-scale covariance matrix estimation and implications for functional
#' genomics.  Statist. Appl. Genet. Mol. Biol. 4:32.
#' @keywords methods multivariate
#' @examples
#' 
#' # Simulate dataset
#' x <- matrix(rnorm(20*30),20,30)
#' # Try different targets
#' shrink.estim(x, tar = build.target(x, type="D"))
#' shrink.estim(x, tar = build.target(x, type="D"))
#' 
#' @importFrom stats cov cor cov2cor
#' @export

shrink.estim <- function(x,tar) {

    if (is.matrix(x)==TRUE && is.numeric(x)==FALSE) stop("The data matrix must be numeric!")

    p <- ncol(x) ; n <- nrow(x)
    covm <- cov(x) ; corm <- cor(x)
    xs <- scale(x,center=TRUE,scale=TRUE)
    
    v <- (n/((n-1)^3))*( crossprod(xs^2) - 1/n*(crossprod(xs))^2 )
    diag(v) <- 0
    
    m <- matrix(rep(apply(xs**2,2,mean),p),p,p)
    f <- (n/(2*(n-1)^3))*(  crossprod(xs**3,xs) + crossprod(xs,xs**3) - (m+t(m))*crossprod(xs)  )
    diag(f) <- 0 ; f[tar == 0] <- 0
    
    corapn <-  cov2cor(tar)
    d      <- (corm - corapn)^2
    lambda <- (sum(v)- sum(corapn*f))/sum(d)
    lambda <- max(min(lambda, 1), 0)
    shrink.cov <-lambda*tar+(1-lambda)*covm
    
    return(list(shrink.cov, c("The shrinkage intensity lambda is:",round(lambda,digits=4))))
}

