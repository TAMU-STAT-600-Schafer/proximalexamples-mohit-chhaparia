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

