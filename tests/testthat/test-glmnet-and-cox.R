test_that("fit.glmnet (binomial) produces a valid glmnet fit", {
  skip_if_not_installed("glmnet")
  set.seed(123)
  n <- 40; p <- 8
  x <- matrix(rnorm(n * p), n, p)
  beta <- c(1, -1, rep(0, p-2))
  eta <- drop(x %*% beta)
  pr <- 1/(1+exp(-eta))
  y <- rbinom(n, 1, pr)

  # Choose a reasonable lambda from an initial path
  g <- glmnet::glmnet(x, y, family = "binomial", standardize = TRUE, nlambda = 20L)
  lam <- g$lambda[round(length(g$lambda) / 2)]

  fit <- fit.glmnet(response = y, x = x, cplx = lam, family = "binomial", standardize = TRUE)
  expect_s3_class(fit, "lognet")

  cf <- coef.glmnet(g, s = 0.1)
  expect_equal(nrow(cf), p + 1)
  expect_true(ncol(cf) >= 1)
  prd <- as.numeric(predict.glmnet(fit, newx = x, type = "link")[, 1])
  expect_length(prd, n)
  expect_true(all(is.finite(prd)))
})

test_that("PLL.coxnet prefers true over permuted outcome", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("survival")
  set.seed(124)
  n <- 60; p <- 6
  x <- matrix(rnorm(n * p), n, p)
  beta <- c(1, 0.8, rep(0, p-2))
  linpred <- drop(scale(x %*% beta))
  u <- runif(n)
  time <- -log(u) / exp(linpred)
  censor <- rexp(n, rate = 1)
  status <- as.integer(time <= censor)
  ytime <- pmin(time, censor)

  fit_cox <- glmnet::glmnet(x, survival::Surv(ytime, status), family = "cox", nlambda = 12L, standardize = TRUE)
  expect_true(inherits(fit_cox, "coxnet"))
  lam <- fit_cox$lambda[round(length(fit_cox$lambda) / 2)]

  pll_true_vec <- c060::PLL.coxnet(fit_cox, newdata = x, newtime = ytime, newstatus = status, complexity = lam)
  pll_perm_vec <- c060::PLL.coxnet(fit_cox, newdata = x, newtime = ytime, newstatus = sample(status), complexity = lam)
  expect_true(is.numeric(pll_true_vec) && is.numeric(pll_perm_vec))
  expect_type(pll_true_vec, "double")
  expect_gt(sum(pll_true_vec), sum(pll_perm_vec) - 1e-6)
})

test_that("tune.glmnet.interval returns a model list with alpha/lambda", {
  skip_if_not_installed("glmnet")
  set.seed(125)
  n <- 50; p <- 10
  x <- matrix(rnorm(n * p), n, p)
  beta <- c(1.2, -0.8, rep(0, p-2))
  eta <- drop(scale(x %*% beta))
  pr <- 1/(1+exp(-eta))
  y <- factor(rbinom(n, 1, pr), levels = c(0,1))
  res <- suppressWarnings(tune.glmnet.interval(parms = 1, x = x, y = y,
                              nfolds = 3L, type.measure = "auc",
                              family = "binomial", verbose = FALSE))
  expect_true(is.list(res))
  expect_true(is.list(res$model))
  expect_true(all(c("alpha", "lambda") %in% names(res$model)))
})


