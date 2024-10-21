# Standardize X and Y
# X - n x p matrix of covariates
# Y - n x 1 response vector
standardizeXY <- function(X, Y){
  n <- nrow(X)
  # Center Y
  Ymean <- mean(Y)
  Ytilde <- Y - Ymean
  # Center X
  Xmeans <- colMeans(X)
  Xprime <- X - matrix(colMeans(X), n, ncol(X), byrow = T)
  # Scale X
  weights <- sqrt(colSums(Xprime^2) / n)
  Xtilde <- Xprime %*% diag(1 / weights)
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}

#  soft-thresholding function
# a - argument to threshold
# lambda - value for threshold
soft <- function(a, lambda){
  return(as.numeric(sign(a) * pmax(0, abs(a) - lambda)))
}

# lasso objective function calculation
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# beta - p vector at which to calculate the function value
# lamdba - tuning parameter
lasso <- function(Xtilde, Ytilde, beta, lambda){
  n <- nrow(Xtilde)
  return(crossprod(Ytilde - Xtilde %*% beta) / (2 * n) + lambda * sum(abs(beta)))
}

# Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta0 - p vector of starting point for coordinate-descent algorithm, optional
# eps - precision level for convergence assessment, default 0.0001
fitLASSOstandardized <- function(Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.0001){
  # Compatibility and other checks
  p <- ncol(Xtilde)
  n <- nrow(Xtilde)
  
  # Coordinate-descent implementation
  error <- 1000
  # Full residual
  r <- Ytilde - Xtilde %*% beta_start
  beta <- beta_start
  while (error > eps){
    beta_old <- beta
    for (j in 1:p){
      # Update of beta
      beta[j] <- soft(beta_old[j] + crossprod(r, Xtilde[,j])/n, lambda)
      # Update of residual
      r <- r - Xtilde[,j] * (beta[j] - beta_old[j])
    }
    # Difference in objective function values
    error <- lasso(Xtilde, Ytilde, beta_old, lambda) - lasso(Xtilde, Ytilde, beta, lambda)
  }
  return(list(beta = beta, fmin = lasso(Xtilde, Ytilde, beta, lambda)))
}



