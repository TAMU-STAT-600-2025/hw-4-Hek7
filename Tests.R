# Creates various test cases for the function to test basic functionality

library(testthat)

source("LassoFunctions.R")

test_that("standardizeXY centers and scales", {
  set.seed(42)
  n <- 100
  p <- 5
  X <- matrix(rnorm(n * p, mean = 5, sd = 2), n, p)
  Y <- matrix(rnorm(n, mean = 10, sd = 3), n, 1)
  
  result <- standardizeXY(X, Y)
  
  # Check center
  expect_equal(mean(result$Ytilde), 0, tolerance = 1e-10)
  expect_equal(result$Ymean, mean(Y))

  expect_equal(colMeans(result$Xtilde), rep(0, p), tolerance = 1e-10)
  
  # Check standardized
  expect_equal(apply(result$Xtilde, 2, sd), rep(1, p), tolerance = 1e-10)
  
  # Check Xmeans
  expect_equal(result$Xmeans, colMeans(X))
  
  # Check weights 
  xcentered <- scale(X, center = TRUE, scale = FALSE)
  expected_weights <- sqrt(colSums(xcentered^2) / n)
  expect_equal(result$weights, expected_weights)
})


test_that("soft thresholding works", {
  # Positive
  expect_equal(soft(5, 2), 3)
  expect_equal(soft(5, 5), 0)
  expect_equal(soft(5, 6), 0)
  
  # Negative
  expect_equal(soft(-5, 2), -3)
  expect_equal(soft(-5, 5), 0)
  expect_equal(soft(-5, 6), 0)
  
  # Zero
  expect_equal(soft(0, 2), 0)
  
  # Vector
  expect_equal(soft(c(5, -3, 1, 0), 2), c(3, -1, 0, 0))
})

test_that("lasso objective function", {
  set.seed(42)
  n <- 50
  p <- 10
  Xtilde <- matrix(rnorm(n * p), n, p)
  Ytilde <- matrix(rnorm(n), n, 1)
  beta <- rnorm(p)
  lambda <- 0.5
  
  result <- lasso(Xtilde, Ytilde, beta, lambda)
  
  # Manual
  residuals <- Ytilde - Xtilde %*% beta
  expected <- sum(residuals^2) / (2 * n) + lambda * sum(abs(beta))
  
  expect_equal(result, expected)
})


test_that("fitLASSOstandardized works with good inputs", {
  set.seed(42)
  n <- 50
  p <- 10
  Xtilde <- matrix(rnorm(n * p), n, p)
  Ytilde <- matrix(rnorm(n), n, 1)
  lambda <- 0.5
  
  result <- fitLASSOstandardized(Xtilde, Ytilde, lambda)
  
  expect_type(result, "list")
  expect_named(result, c("beta", "fmin"))
  expect_length(result$beta, p)
  expect_length(result$fmin, 1)
  expect_true(result$fmin >= 0)
})

test_that("fitLASSOstandardized catches bad inputs", {
  set.seed(42)
  n <- 50
  p <- 10
  Xtilde <- matrix(rnorm(n * p), n, p)
  Ytilde <- matrix(rnorm(n), n, 1)
  
  # dim error
  Ytilde_wrong <- matrix(rnorm(n + 5), n + 5, 1)
  expect_error(fitLASSOstandardized(Xtilde, Ytilde_wrong, 0.5))
  
  # <0 lambda
  expect_error(fitLASSOstandardized(Xtilde, Ytilde, -0.5))
  
  # bad beta dim
  beta_wrong <- rnorm(p + 2)
  expect_error(fitLASSOstandardized(Xtilde, Ytilde, 0.5, beta_start = beta_wrong))
})

test_that("fitLASSOstandardized uses starting values", {
  set.seed(123)
  n <- 50
  p <- 10
  Xtilde <- matrix(rnorm(n * p), n, p)
  Ytilde <- matrix(rnorm(n), n, 1)
  lambda <- 0.5
  
  beta_start <- rnorm(p)
  result <- fitLASSOstandardized(Xtilde, Ytilde, lambda, beta_start = beta_start)
  
  expect_length(result$beta, p)
})

test_that("fitLASSOstandardized gives zero solution lambda >> 0", {
  set.seed(123)
  n <- 50
  p <- 10
  Xtilde <- matrix(rnorm(n * p), n, p)
  Ytilde <- matrix(rnorm(n), n, 1)
  
  lambda_max <- max(abs(crossprod(Xtilde, Ytilde))) / n
  result <- fitLASSOstandardized(Xtilde, Ytilde, lambda_max * 2)
  
  expect_equal(result$beta, rep(0, p), tolerance = 1e-6)
})

# Test fitLASSOstandardized_seq
test_that("fitLASSOstandardized_seq works correctly", {
  set.seed(123)
  n <- 50
  p <- 10
  Xtilde <- matrix(rnorm(n * p), n, p)
  Ytilde <- matrix(rnorm(n), n, 1)
  n_lambda <- 20
  
  result <- fitLASSOstandardized_seq(Xtilde, Ytilde, n_lambda = n_lambda)
  
  expect_type(result, "list")
  expect_named(result, c("lambda_seq", "beta_mat", "fmin_vec"))
  expect_length(result$lambda_seq, n_lambda)
  expect_equal(dim(result$beta_mat), c(p, n_lambda))
  expect_length(result$fmin_vec, n_lambda)
  
  # Lambda grows smaller 
  expect_true(all(diff(result$lambda_seq) < 0))
})

test_that("fitLASSOstandardized_seq validates", {
  set.seed(123)
  n <- 50
  p <- 10
  Xtilde <- matrix(rnorm(n * p), n, p)
  Ytilde <- matrix(rnorm(n), n, 1)
  
  # bad dims
  Ytilde_wrong <- matrix(rnorm(n + 5), n + 5, 1)
  expect_error(fitLASSOstandardized_seq(Xtilde, Ytilde_wrong))
})

test_that("fitLASSOstandardized_seq handles custom lambda sequences", {
  set.seed(123)
  n <- 50
  p <- 10
  Xtilde <- matrix(rnorm(n * p), n, p)
  Ytilde <- matrix(rnorm(n), n, 1)
  
  lambda_seq <- c(1.0, 0.5, 0.1, 0.01)
  result <- fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq = lambda_seq)
  
  expect_equal(result$lambda_seq, lambda_seq)
  expect_equal(ncol(result$beta_mat), length(lambda_seq))
})

test_that("fitLASSOstandardized_seq filters lambda < 0  ", {
  set.seed(123)
  n <- 50
  p <- 10
  Xtilde <- matrix(rnorm(n * p), n, p)
  Ytilde <- matrix(rnorm(n), n, 1)
  
  lambda_seq <- c(1.0, 0.5, -0.1, 0.01)
  result <- fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq = lambda_seq)
  
  expect_true(all(result$lambda_seq >= 0))
  expect_equal(result$lambda_seq, c(1.0, 0.5, 0.01))
})

test_that("fitLASSOstandardized_seq warns on all negative lambdas", {
  set.seed(123)
  n <- 50
  p <- 10
  Xtilde <- matrix(rnorm(n * p), n, p)
  Ytilde <- matrix(rnorm(n), n, 1)
  
  lambda_seq <- c(-1.0, -0.5, -0.1)
  expect_warning(fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq = lambda_seq))
})

test_that("fitLASSO works on original data", {
  set.seed(123)
  n <- 50
  p <- 10
  X <- matrix(rnorm(n * p, mean = 5), n, p)
  Y <- matrix(rnorm(n, mean = 10), n, 1)
  
  result <- fitLASSO(X, Y, n_lambda = 20)
  
  expect_type(result, "list")
  expect_named(result, c("lambda_seq", "beta_mat", "beta0_vec"))
  expect_equal(nrow(result$beta_mat), p)
  expect_length(result$beta0_vec, length(result$lambda_seq))
})

test_that("fitLASSO back-transforms coefficients correctly", {
  set.seed(123)
  n <- 100
  p <- 5
  X <- matrix(rnorm(n * p, mean = 5, sd = 2), n, p)
  beta_true <- c(1, -1, 0.5, 0, 0)
  Y <- matrix(X %*% beta_true + rnorm(n, sd = 0.1), n, 1)
  
  result <- fitLASSO(X, Y, n_lambda = 50)
  
  lambda_idx <- length(result$lambda_seq) 
  Y_pred <- result$beta0_vec[lambda_idx] + X %*% result$beta_mat[, lambda_idx]
  
  expect_true(cor(Y, Y_pred) > 0.9)
})

test_that("cvLASSO performs cross-validation", {
  set.seed(123)
  n <- 100
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n), n, 1)
  
  result <- cvLASSO(X, Y, n_lambda = 20, k = 5)
  
  expect_type(result, "list")
  expect_named(result, c("lambda_seq", "beta_mat", "beta0_vec", "fold_ids", 
                         "lambda_min", "lambda_1se", "cvm", "cvse"))
  
  expect_length(result$fold_ids, n)
  expect_true(all(result$fold_ids %in% 1:5))
  
  expect_true(result$lambda_min %in% result$lambda_seq)
  expect_true(result$lambda_1se %in% result$lambda_seq)
  expect_true(result$lambda_1se >= result$lambda_min)
  
  expect_length(result$cvm, length(result$lambda_seq))
  expect_length(result$cvse, length(result$lambda_seq))
  expect_true(all(result$cvm >= 0))
  expect_true(all(result$cvse >= 0))
})

test_that("cvLASSO uses custom fold_ids", {
  set.seed(123)
  n <- 100
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n), n, 1)
  
  custom_folds <- sample(rep(1:3, length.out = n))
  result <- cvLASSO(X, Y, fold_ids = custom_folds)
  
  expect_equal(result$fold_ids, custom_folds)
  expect_equal(max(result$fold_ids), 3)
})

test_that("cvLASSO valid Lambda", {
  set.seed(123)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- c(rep(1, 5), rep(0, 15))
  Y <- matrix(X %*% beta_true + rnorm(n, sd = 0.5), n, 1)
  
  result <- cvLASSO(X, Y, n_lambda = 50, k = 5)
  
  min_idx <- which(result$lambda_seq == result$lambda_min)
  max_idx <- 1
  expect_true(result$cvm[min_idx] <= result$cvm[max_idx])
  
  expect_true(result$lambda_1se >= result$lambda_min)
})

test_that("All Functions Together", {
  set.seed(123)
  n <- 200
  p <- 30
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- c(rep(2, 3), rep(0, 27))
  Y <- matrix(X %*% beta_true + rnorm(n), n, 1)
  
  cv_result <- cvLASSO(X, Y, n_lambda = 40, k = 5)
  
  lambda_min_idx <- which(cv_result$lambda_seq == cv_result$lambda_min)
  beta_hat <- cv_result$beta_mat[, lambda_min_idx]
  beta0_hat <- cv_result$beta0_vec[lambda_min_idx]
  
  Y_pred <- beta0_hat + X %*% beta_hat
  
  expect_true(cor(Y, Y_pred) > 0.5)
  
  expect_true(sum(abs(beta_hat) > 0.1) >= 1)
})

