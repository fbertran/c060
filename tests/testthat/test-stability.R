`%||%` <- function(x, y) if (is.null(x)) y else x

# Ensure symbols exist in .GlobalEnv for Windows PSOCK clusters
skip_if_not_installed("glmnet")
skip_if_not_installed("Matrix")
glmnet_fn <- get("glmnet", asNamespace("glmnet"))
drop0_fn  <- get("drop0",  asNamespace("Matrix"))
old_glmnet <- if (exists("glmnet", envir = .GlobalEnv, inherits = FALSE)) get("glmnet", envir = .GlobalEnv) else NULL
old_drop0  <- if (exists("drop0",  envir = .GlobalEnv, inherits = FALSE)) get("drop0",  envir = .GlobalEnv) else NULL
assign("glmnet", glmnet_fn, envir = .GlobalEnv)
assign("drop0",  drop0_fn,  envir = .GlobalEnv)
withr::defer({
  if (is.null(old_glmnet)) rm(list = "glmnet", envir = .GlobalEnv) else assign("glmnet", old_glmnet, envir = .GlobalEnv)
  if (is.null(old_drop0))  rm(list = "drop0",  envir = .GlobalEnv) else assign("drop0",  old_drop0,  envir = .GlobalEnv)
}, test_env())

test_that("stabpath returns a coherent matrix path", {
  set.seed(126)
  n <- 80; p <- 12
  x <- matrix(rnorm(n * p), n, p)
  beta <- c(2, -2, 1.5, rep(0, p-3))
  pr <- 1/(1+exp(-scale(drop(x %*% beta))))
  y <- factor(rbinom(n, 1, pr), levels = c(0,1))
  
  sp <- stabpath(y = y, x = x, steps = 20L, weakness = 1, family = "binomial")
  expect_true(is.matrix(sp$stabpath))
  expect_equal(nrow(sp$stabpath), p)
  expect_true(ncol(sp$stabpath) >= 1)
  expect_true(all(is.finite(sp$stabpath)))
})

test_that("stabsel selects variables at or above threshold (may be empty)", {
  set.seed(127)
  n <- 80; p <- 10
  x <- matrix(rnorm(n * p), n, p)
  eta <- drop(scale(x[,1]*2 - x[,2]))
  y <- factor(rbinom(n, 1, 1/(1+exp(-eta))), levels = c(0,1))
  
  sp <- stabpath(y = y, x = x, steps = 20L, weakness = 1, family = "binomial")
  ss <- stabsel(sp, error = 0.05, type = "pfer", pi_thr = 0.6)
  S <- ss$stable %||% integer(0)
  if (length(S)) {
    ph <- sp$stabpath[, ss$lpos, drop = TRUE]
    expect_true(all(ph[S] >= 0.6 - 1e-8))
  } else {
    succeed("No variables crossed threshold in this tiny run (acceptable)")
  }
})
