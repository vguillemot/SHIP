test_that("all target builders return finite covariance matrices", {
  set.seed(1)
  x <- matrix(rnorm(40), nrow = 10, ncol = 4)
  groups <- list("a", "a", "b", character())

  for (type in c("D", "F", "G", "Gpos", "Gstar", "cor")) {
    target <- build.target(x, groups, type)
    expect_equal(dim(target), c(4L, 4L))
    expect_true(all(is.finite(target)))
    expect_equal(target, t(target))
  }
})

test_that("targetCor handles a single variable", {
  set.seed(3)
  x <- matrix(rnorm(10), ncol = 1L)
  expect_equal(targetCor(x, list("a")), matrix(var(x), 1L, 1L))
})

test_that("targetCor averages only significant linked correlations", {
  x <- cbind(seq_len(10), seq_len(10), rep(c(-1, 1), 5))
  groups <- list("a", "a", "a")
  target <- targetCor(x, groups)

  expect_equal(target[1, 2], var(x[, 1]))
  expect_equal(target[1, 3], 0)
  expect_equal(target[2, 3], 0)
})

test_that("target helpers handle empty groups and invalid inputs", {
  expect_equal(target.help(list()), matrix(numeric(), 0L, 0L))
  expect_equal(check.path("path", "path"), 1L)
  expect_equal(check.path("path", "other"), 0L)
  expect_error(target.help(NULL), "must be a list")
  expect_error(target.help(list(NULL, "path")), "each genegroups entry")
  expect_error(target.help(list(list("path"), "path")), "each genegroups entry")

  x <- matrix(rnorm(20), nrow = 10, ncol = 2)
  expect_error(build.target(x, type = "unknown"))
  expect_error(build.target(x, type = "G"), "genegroups")
  expect_error(build.target(matrix(1, 2, 1), type = "D"), "non-zero variance")
})

test_that("shrink.estim returns a named numeric intensity", {
  set.seed(2)
  x <- matrix(rnorm(30), nrow = 10, ncol = 3)
  fit <- shrink.estim(x, cov(x))

  expect_named(fit, c("shrink.cov", "lambda"))
  expect_true(is.numeric(fit$lambda))
  expect_equal(fit$lambda, 0)
  expect_equal(fit$shrink.cov, cov(x))
})
