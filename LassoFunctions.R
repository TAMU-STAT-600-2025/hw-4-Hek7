# [ToDo] Standardize X and Y: center both X and Y; scale centered X
# X - n x p matrix of covariates
# Y - n x 1 response vector
standardizeXY <- function(X, Y){
  # [ToDo] Center Y
  Ymean <- mean(Y)
  Ytilde <- Y - y_mean
  
  # [ToDo] Center and scale X
  Xmeans <- colMeans(X)
  
  xcentered <- scale(X, center = TRUE, scale = FALSE) 
  Xtilde <- scale(X, center = TRUE, scale = TRUE) 
  
  n <- nrow(X)
  weights <- sqrt(colSums(xcentered^2) / n)
  
  # Return:
  # Xtilde - centered and appropriately scaled X
  # Ytilde - centered Y
  # Ymean - the mean of original Y
  # Xmeans - means of columns of X (vector)
  # weights - defined as sqrt(X_j^{\top}X_j/n) after centering of X but before scaling
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}

# [ToDo] Soft-thresholding of a scalar a at level lambda 
# [OK to have vector version as long as works correctly on scalar; will only test on scalars]
soft <- function(a, lambda){
  return(sign(a)*pmax(abs(a) - lambda, 0))
}

# [ToDo] Calculate objective function of lasso given current values of Xtilde, Ytilde, beta and lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta - value of beta at which to evaluate the function
lasso <- function(Xtilde, Ytilde, beta, lambda){
  n = nrow(Xtilde)
  
  # Calculate objective
  fobj = ((2*n)^(-1)) * norm(Ytilde - Xtilde%*%beta, type = c("2"))^2 + lambda*sum(abs(beta))
}

# [ToDo] Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1 (vector)
# lamdba - tuning parameter
# beta_start - p vector, an optional starting point for coordinate-descent algorithm
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized <- function(Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.001){
  #[ToDo]  Check that n is the same between Xtilde and Ytilde
  
  
  if (nrow(Xtilde) != nrow(Ytilde)){
    stop("rows in ytilde and xtilde not equal")
  }
  
  n = nrow(Xtilde)
  p = ncol(Xtilde)
  
  #[ToDo]  Check that lambda is non-negative
  if(lambda < 0){
    stop("lambda must be non-negative")
  }
  
  
  #[ToDo]  Check for starting point beta_start. 
  # If none supplied, initialize with a vector of zeros.
  # If supplied, check for compatibility with Xtilde in terms of p
  
  if(is.null(beta_start)){
    beta <- rep(0, p)
  } else {
    if(length(beta_start) != p){
      stop("beta_start dim does not match number of columns in Xtilde")
    }
    beta <- beta_start
  }
  
  #[ToDo]  Coordinate-descent implementation. 
  # Stop when the difference between objective functions is less than eps for the first time.
  # For example, if you have 3 iterations with objectives 3, 1, 0.99999,
  # your should return fmin = 0.99999, and not have another iteration
  
  fobj_old <- lasso(Xtilde, Ytilde, beta, lambda)
  converged <- FALSE
  
  while(!converged){
    for(j in 1:p){
      r_j <- Ytilde - Xtilde[, -j, drop = FALSE] %*% beta[-j]
      beta[j] <- soft((1/n) * crossprod(Xtilde[, j], r_j), lambda)
    }
    
    # new obj
    fobj_new <- lasso(Xtilde, Ytilde, beta, lambda)
    
    # Check condition
    if(abs(fobj_old - fobj_new) < eps){
      converged <- TRUE
    }
    
    # Update 
    fobj_old <- fobj_new
  }
  
  # save best 
  fmin <- fobj_new
  
  # Return 
  # beta - the solution (a vector)
  # fmin - optimal function value (value of objective at beta, scalar)
  return(list(beta = beta, fmin = fmin))
}

# [ToDo] Fit LASSO on standardized data for a sequence of lambda values. Sequential version of a previous function.
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence,
#             is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized_seq <- function(Xtilde, Ytilde, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  # [ToDo] Check that n is the same between Xtilde and Ytilde
  if (nrow(Xtilde) != nrow(Ytilde)){
    stop("rows in ytilde and xtilde not equal")
  }
  
  n = nrow(Xtilde)
  p = ncol(Xtilde)
 
  # [ToDo] Check for the user-supplied lambda-seq (see below)
  # If lambda_seq is supplied, only keep values that are >= 0,
  # and make sure the values are sorted from largest to smallest.
  # If none of the supplied values satisfy the requirement,
  # print the warning message and proceed as if the values were not supplied.
  
  if(!is.null(lambda_seq)){
    lambda_seq <- lambda_seq[lambda_seq >= 0]
    
    # Check if any valid values remain
    if(length(lambda_seq) == 0){
      warning("No non-negative lambda values supplied. Generating default sequence.")
      lambda_seq <- NULL
    } else {
      lambda_seq <- sort(lambda_seq, decreasing = TRUE)
    }
  }
  
  
  # If lambda_seq is not supplied, calculate lambda_max 
  # (the minimal value of lambda that gives zero solution),
  # and create a sequence of length n_lambda as
  
  if(is.null(lambda_seq)){
    lambda_max <- max(abs(crossprod(Xtilde, Ytilde))) / n
    
    lambda_seq <- exp(seq(log(lambda_max), log(0.01), length = n_lambda))
  }
  
  # [ToDo] Apply fitLASSOstandardized going from largest to smallest lambda 
  # (make sure supplied eps is carried over). 
  # Use warm starts strategy discussed in class for setting the starting values.
  
  n_seq <- length(lambda_seq)
  beta_mat <- matrix(0, nrow = p, ncol = n_seq)
  fmin_vec <- numeric(n_seq)
  
  beta_start <- NULL # defaults to 0's in the fit func
  
  for(i in 1:n_seq){
    fit <- fitLASSOstandardized(Xtilde, Ytilde, lambda_seq[i], 
                                beta_start = beta_start, eps = eps)

    beta_mat[, i] <- fit$beta
    fmin_vec[i] <- fit$fmin
    
    beta_start <- fit$beta
  }
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value
  # fmin_vec - length(lambda_seq) vector of corresponding objective function values at solution
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, fmin_vec = fmin_vec))
}

# [ToDo] Fit LASSO on original data using a sequence of lambda values
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  # [ToDo] Center and standardize X,Y based on standardizeXY function
  standardized <- standardizeXY(X, Y)
  Xtilde <- standardized$Xtilde
  Ytilde <- standardized$Ytilde
  Ymean <- standardized$Ymean
  Xmeans <- standardized$Xmeans
  weights <- standardized$weights
  
  # [ToDo] Fit Lasso on a sequence of values using fitLASSOstandardized_seq
  # (make sure the parameters carry over)
  fit <- fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq = lambda_seq, 
                                          n_lambda = n_lambda, eps = eps)
  
  lambda_seq <- fit$lambda_seq
  beta_mat_std <- fit$beta_mat
 
  # [ToDo] Perform back scaling and centering to get original intercept and coefficient vector
  # for each lambda
  p <- ncol(X)
  n_lambdas <- length(lambda_seq)

  beta_mat <- matrix(0, nrow = p, ncol = n_lambdas)
  beta0_vec <- numeric(n_lambdas)
  
  for(i in 1:n_lambdas){
    beta_mat[, i] <- beta_mat_std[, i] / weights
    
    beta0_vec[i] <- Ymean - sum(Xmeans * beta_mat[, i])
  }
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec))
}


# [ToDo] Fit LASSO and perform cross-validation to select the best fit
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# k - number of folds for k-fold cross-validation, default is 5
# fold_ids - (optional) vector of length n specifying the folds assignment (from 1 to max(folds_ids)), if supplied the value of k is ignored 
# eps - precision level for convergence assessment, default 0.001
cvLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, k = 5, fold_ids = NULL, eps = 0.001){
  n <- nrow(X)
  
  # [ToDo] Fit Lasso on original data using fitLASSO
  fit_orig <- fitLASSO(X, Y, lambda_seq = lambda_seq, n_lambda = n_lambda, eps = eps)
  
  lambda_seq <- fit_orig$lambda_seq
  n_lambdas <- length(lambda_seq)
  beta_mat <- fit_orig$beta_mat
  beta0_vec <- fit_orig$beta0_vec
 
  # [ToDo] If fold_ids is NULL, split the data randomly into k folds.
  # If fold_ids is not NULL, split the data according to supplied fold_ids.
  
  if(is.null(fold_ids)){
    fold_ids <- sample(rep(1:k, length.out = n))
  } else {
    k <- max(fold_ids)
  }
  
  # [ToDo] Calculate LASSO on each fold using fitLASSO,
  # and perform any additional calculations needed for CV(lambda) and SE_CV(lambda)
  mse_matrix <- matrix(0, nrow = k, ncol = n_lambdas)
  
  for(fold in 1:k){
    val_idx <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)
    
    X_train <- X[train_idx, , drop = FALSE]
    Y_train <- Y[train_idx, , drop = FALSE]
    X_val <- X[val_idx, , drop = FALSE]
    Y_val <- Y[val_idx, , drop = FALSE]
    
    fold_fit <- fitLASSO(X_train, Y_train, lambda_seq = lambda_seq, eps = eps)
    
    for(i in 1:n_lambdas){
      Y_pred <- fold_fit$beta0_vec[i] + X_val %*% fold_fit$beta_mat[, i]
      mse_matrix[fold, i] <- mean((Y_val - Y_pred)^2)
    }
  }
  
  # Calculate CV(lambda) - mean MSE across folds for each lambda
  cvm <- colMeans(mse_matrix)
  
  # Calculate SE_CV(lambda) - standard error of MSE across folds
  cvse <- apply(mse_matrix, 2, function(x) sd(x) / sqrt(k))
  
  # [ToDo] Find lambda_min
  min_idx <- which.min(cvm)
  lambda_min <- lambda_seq[min_idx]
  
  # [ToDo] Find lambda_1SE
  min_cvm <- cvm[min_idx]
  min_cvse <- cvse[min_idx]
  threshold <- min_cvm + min_cvse
  
  lambda_1se_idx <- min(which(cvm <= threshold))
  lambda_1se <- lambda_seq[lambda_1se_idx]
  
  # Return output
  # Output from fitLASSO on the whole data
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  # fold_ids - used splitting into folds from 1 to k (either as supplied or as generated in the beginning)
  # lambda_min - selected lambda based on minimal rule
  # lambda_1se - selected lambda based on 1SE rule
  # cvm - values of CV(lambda) for each lambda
  # cvse - values of SE_CV(lambda) for each lambda
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec, fold_ids = fold_ids, lambda_min = lambda_min, lambda_1se = lambda_1se, cvm = cvm, cvse = cvse))
}

