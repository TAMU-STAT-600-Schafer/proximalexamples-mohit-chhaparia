# Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta0 - p vector of starting point for coordinate-descent algorithm, optional
# eps - precision level for convergence assessment, default 0.0001
# (NEW) s - step size for proximal gradient
fitLASSOstandardized_prox <- function(Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.0001, s = 0.01){
  # Compatibility and other checks
  p <- ncol(Xtilde)
  n <- nrow(Xtilde)
  if (length(Ytilde) != n){
    stop("Dimensions of X and Y don't match")
  }
  if (lambda < 0){
    stop("Only non-negative values of lambda are allowed!")
  }
  if (is.null(beta_start)){
    # Initialize beta0
    beta_start <- rep(0,p)
  }else if (length(beta_start) != p){
    stop("Supplied initial starting point has length different from p!")
  }
  
  # Proximal gradient-descent implementation
  
  error <- 1000
  
  # Calculate full residual
  r <- Ytilde - Xtilde %*% beta_start
  beta <- beta_start
  
  while(error > eps){
    
    beta_old <- beta
    
    # move in direction of gradient of MSE
    beta_move <- beta + s * crossprod(Xtilde, r) / n
    
    # apply prox operator (soft) to beta_move
    beta <- soft(beta_move, s * lambda)
    
    # update the residual
    r <- r - Xtilde %*% (beta - beta_old)
    
    # calculate error
    error <- abs(lasso(Xtilde, Ytilde, beta_old, lambda) - 
      lasso(Xtilde, Ytilde, beta, lambda))
  }
  
  return(list(beta = beta, fmin = lasso(Xtilde, Ytilde, beta, lambda)))
}



# Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta0 - p vector of starting point for coordinate-descent algorithm, optional
# eps - precision level for convergence assessment, default 0.0001
# (NEW) s - step size for proximal gradient
fitLASSOstandardized_prox_Nesterov <- function(Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.0001, s = 0.01){
  # Compatibility and other checks
  p <- ncol(Xtilde)
  n <- nrow(Xtilde)
  if (length(Ytilde) != n){
    stop("Dimensions of X and Y don't match")
  }
  if (lambda < 0){
    stop("Only non-negative values of lambda are allowed!")
  }
  if (is.null(beta_start)){
    # Initialize beta0
    beta_start <- rep(0,p)
  }else if (length(beta_start) != p){
    stop("Supplied initial starting point has length different from p!")
  }
  
  # Proximal gradient-descent implementation
  
  error <- 1000
  
  # Calculate full residual
  r <- Ytilde - Xtilde %*% beta_start
  beta <- beta_start
  beta <- x_old <- beta_start
  
  # Initialize the Nesterov Sequence
  
  lambda_t <- 1
  lambda_t1 <- (1 + sqrt(1 + 4 * lambda_t ^ 2)) / 2
  
  while(error > eps){
    
    beta_old <- beta
    
    # move in direction of gradient of MSE
    beta_move <- beta + s * crossprod(Xtilde, r) / n
    
    # apply prox operator (soft) to beta_move
    x <- soft(beta_move, s * lambda)
    
    # apply the Nesterov acceleration
    beta <- x + ((lambda_t - 1) / lambda_t1) * (x - x_old)
    
    # update the residual
    r <- r - Xtilde %*% (beta - beta_old)
    
    # update x_old
    x_old <- x
    
    #update lambda_t and lambda_t1
    lambda_t <- lambda_t1
    lambda_t1 <- (1 + sqrt(1 + 4 * lambda_t ^ 2)) / 2
    
    # calculate error
    error <- abs(lasso(Xtilde, Ytilde, beta_old, lambda) - 
                   lasso(Xtilde, Ytilde, beta, lambda))
  }
  
  return(list(beta = beta, fmin = lasso(Xtilde, Ytilde, beta, lambda)))
}


# Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lambda - tuning parameter
# beta0 - p vector of starting point for coordinate-descent algorithm, optional
# eps - precision level for convergence assessment, default 0.0001
# (NEW) tau - ADMM parameter
fitLASSOstandardized_ADMM <- function(Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.0001, tau = 0.01){
  # Compatibility and other checks
  p <- ncol(Xtilde)
  n <- nrow(Xtilde)
  if (length(Ytilde) != n){
    stop("Dimensions of X and Y don't match")
  }
  if (lambda < 0){
    stop("Only non-negative values of lambda are allowed!")
  }
  if (is.null(beta_start)){
    # Initialize beta0
    beta_start <- rep(0,p)
  }else if (length(beta_start) != p){
    stop("Supplied initial starting point has length different from p!")
  }
  
  # ADMM implementation


  return(list(beta = beta, fmin = lasso(Xtilde, Ytilde, beta, lambda)))
}

