#' SHIP provides shrinkage covariance estimation with user-selected targets.
#' The available targets include diagonal, constant-correlation, and
#' knowledge-based structures informed by functional gene groups.
#'
#' Start with \code{\link{build.target}} to construct a target matrix, then
#' pass it to \code{\link{shrink.estim}} together with the data matrix.
#' 
#' @author Monika Jelizarow and Vincent Guillemot
#' @references \itemize{ \item J. Schaefer and K. Strimmer, 2005. A shrinkage
#' approach to large-scale covariance matrix estimation and implications for
#' functional genomics.  Statist. Appl. Genet. Mol. Biol. 4:32. \item M.
#' Jelizarow, V. Guillemot, A. Tenenhaus, K. Strimmer, A.-L. Boulesteix, 2010.
#' Over-optimism in bioinformatics: an illustration. Bioinformatics. Accepted.
#' }
#' @examples
#' data("expl")
#' target <- build.target(expl$x, expl$genegroups, type = "G")
#' estimate <- shrink.estim(expl$x, target)
#' estimate$lambda
#' @keywords package
"_PACKAGE"

