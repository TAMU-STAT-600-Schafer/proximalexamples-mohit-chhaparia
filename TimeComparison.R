#########################################
# Source the functions and necessary libraries
#########################################
# Lasso functions for coordinate descent
source("LASSO_CoordinateDescent.R")
# Functions written in the starter code
source("StarterCode.R")
# For later time comparisons
library(microbenchmark)
#########################################
# Generate toy example dataset, standardize
#########################################
set.seed(38947)
p = 50
n = 30
Y <- rnorm(n)
X <- matrix(rnorm(n*p), n, p)
out <- standardizeXY(X,Y)
lambda_max <- max(abs(crossprod(out$Xtilde, out$Ytilde))/nrow(X))

#########################################
# Check equality of outputs between coordinate descent and proximal gradient algorithms
##########################################
lambda1 = lambda_max/2
out_coord <- fitLASSOstandardized(out$Xtilde, out$Ytilde, beta_start = rep(0, p), lambda = lambda1, eps = 1e-10)
out_prox <- fitLASSOstandardized_prox(out$Xtilde, out$Ytilde, beta_start = rep(0, p), lambda = lambda1, eps = 1e-10, s = 0.1)
out_coord$fmin - out_prox$fmin
plot(out_coord$beta, out_prox$beta)


# Check the implementation time
microbenchmark(
  fitLASSOstandardized(out$Xtilde, out$Ytilde, beta_start = rep(0, p),lambda = lambda1, eps = 1e-6),
  fitLASSOstandardized_prox(out$Xtilde, out$Ytilde, beta_start = rep(0, p),lambda = lambda1, eps = 1e-6, s = 0.1)
)

#########################################
# Test Nesterov acceleration
##########################################
out_prox2 <- fitLASSOstandardized_prox_Nesterov(out$Xtilde, out$Ytilde, beta_start = rep(0, p),lambda = lambda1, eps = 1e-10, s = 0.1)
out_coord$fmin - out_prox2$fmin
plot(out_coord$beta, out_prox2$beta)
plot(out_prox2$fobj_vec[-c(1:40)])

# Check the implementation time
microbenchmark(
  fitLASSOstandardized(out$Xtilde, out$Ytilde, beta_start = rep(0, p),lambda = lambda1, eps = 1e-10),
  fitLASSOstandardized_prox(out$Xtilde, out$Ytilde, beta_start = rep(0, p),lambda = lambda1, eps = 1e-10, s = 0.1),
  fitLASSOstandardized_prox_Nesterov(out$Xtilde, out$Ytilde, beta_start = rep(0, p),lambda = lambda1, eps = 1e-10, s = 0.1)
)

#########################################
# Test ADMM implementation
##########################################

out_admm<- fitLASSOstandardized_ADMM(out$Xtilde, out$Ytilde, beta_start = rep(0, p),lambda = lambda1, eps = 1e-10, tau = 0.1)
out_coord$fmin - out_admm$fmin
plot(out_coord$beta, out_admm$beta)

# Check the implementation time
microbenchmark(
  fitLASSOstandardized(out$Xtilde, out$Ytilde, beta_start = rep(0, p),lambda = lambda1, eps = 1e-10),
  fitLASSOstandardized_prox(out$Xtilde, out$Ytilde, beta_start = rep(0, p),lambda = lambda1, eps = 1e-10, s = 0.1),
  fitLASSOstandardized_prox_Nesterov(out$Xtilde, out$Ytilde, beta_start = rep(0, p), lambda = lambda1, eps = 1e-10, s = 0.1),
  fitLASSOstandardized_ADMM(out$Xtilde, out$Ytilde, beta_start = rep(0, p),lambda = lambda1, eps = 1e-10, tau = 0.1)
)
