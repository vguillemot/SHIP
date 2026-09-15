validate_data <- function(x) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("x must be a numeric matrix.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("x must contain only finite values.", call. = FALSE)
  }
  if (nrow(x) < 2L || ncol(x) < 1L) {
    stop("x must have at least two observations and one variable.", call. = FALSE)
  }
  if (any(apply(x, 2L, var) == 0)) {
    stop("x must have non-zero variance in every column.", call. = FALSE)
  }
  invisible(x)
}

is_valid_gene_group <- function(gene_group) {
  !is.null(gene_group) && is.atomic(gene_group) && !is.object(gene_group)
}

validate_gene_groups <- function(gene_groups, n_variables) {
  if (!is.list(gene_groups) || length(gene_groups) != n_variables) {
    stop("genegroups must be a list with one entry per column of x.", call. = FALSE)
  }
  if (any(!vapply(gene_groups, is_valid_gene_group, logical(1)))) {
    stop("each genegroups entry must be an atomic vector of pathway IDs.", call. = FALSE)
  }
  invisible(gene_groups)
}

validate_target <- function(target, n_variables) {
  if (!is.matrix(target) || !is.numeric(target) ||
      !identical(dim(target), c(n_variables, n_variables)) ||
      any(!is.finite(target))) {
    stop("target must be a finite numeric square matrix matching x.", call. = FALSE)
  }
  invisible(target)
}
