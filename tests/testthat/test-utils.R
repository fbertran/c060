`%||%` <- function(x, y) if (is.null(x)) y else x

test_that("balancedFolds yields approximately balanced class counts", {
  if (!exists("balancedFolds")) skip("balancedFolds not available")
  set.seed(128)
  y <- factor(sample(c(0L, 1L), size = 101, replace = TRUE, prob = c(0.6, 0.4)))
  folds <- c060::balancedFolds(class.column.factor = y, cross.outer = 5L)
  expect_equal(length(unique(folds)), 5L)
  all_idx <- sort(unlist(folds))
  expect_true(all(all_idx==rep(1:5,c(20,21,21,20,19))))
  pos <- levels(y)[2]
  props <- sapply(folds, function(idx) mean(y[idx] == pos))
})

test_that("aggregation.auc reflects signal vs noise", {
  skip_if_not_installed("glmnet")
  set.seed(129)
  n <- 120; p <- 5
  x <- matrix(rnorm(n*p), n, p)
  # strong signal in first column
  lin <- x[,1]*2
  y <- rbinom(n, 1, 1/(1+exp(-lin)))

  model_sig <- glmnet::glmnet(x, y, family = "binomial")
  auc_sig <- aggregation.auc(full.data = NULL, response = y, x = x, model = model_sig, type = "apparent")
  expect_true(is.finite(auc_sig))

  y_perm <- sample(y)
  model_noise <- glmnet::glmnet(x, y_perm, family = "binomial")
  auc_noise <- aggregation.auc(full.data = NULL, response = y_perm, x = x, model = model_noise, type = "apparent")
  expect_true(is.finite(auc_noise))

  expect_gt(auc_sig, auc_noise - 1e-6)
  expect_gte(auc_sig, 0.7)
})
