### NEWTON RAPHSON ALGORITHM ###
# Goal: find the root (x axis intersection) of a function, works for real and complex functions
# 1. pick a random point x of the function ("first guess")
# 2. find the tangent to the curve at that point x
# 3. find the intersection of the tangent with the x axis and make it your new x
# 4. repeat until you're close enough (delta-y below a certain threshold which defines that you're on a point)

## example on weibull from litterature
library(ggplot2)
GP_data <-
  data.frame("lifetimes" = c(10, 33, 44, 56, 59, 72, 74, 77, 92, 93, 96, 100, 100, 102, 105, 10))
NR_fit_Weibull <- function(x, beta_0, eps = 0.000001)
{
  beta <- beta_0
  n <- length(x)
  
  sumlogx <- sum(log(x))
  
  diff <- 1
  alpha <- beta_0[1]
  theta <- beta_0[2]
  
  
  while (diff > eps)
  {
    beta.old <- beta
    
    w1 <- sum(x ^ alpha)
    w2 <- sum(x ^ alpha * log(x))
    w3 <- sum(x ^ alpha * (log(x) ^ 2))
    s <- c(n / alpha + sumlogx - w2 / theta,-n /
            theta + w1 / theta ^ 2)
    J <- matrix(c(
      -n / alpha ^ 2 - w3 / theta,
      w2 / theta ^ 2,
      w2 / theta ^ 2,
      n / theta ^ 2 - w1 * 2 / (theta ^ 3)
    ),
    ncol <- 2)
    beta <- beta - solve(J) %*% s
    alpha <- beta[1]
    theta <- beta[2]
    diff <- sum(abs(beta - beta.old))
  }
  list(beta = beta, J = J)
}
lifetimes <- GP_data$lifetimes[GP_data$lifetimes != 0]
NR_fit_Weibull(x = lifetimes,
               beta_0 = c(5, 155),
               eps = 0.001)


dloglogis <- function(x, alpha, beta) {
  res <- (beta/alpha) * (x/alpha)^(beta-1)*(1+(x/alpha)^beta)^-2
  return(res)
}

## for log-logistic ça donne ça
NR_fit_LogLogistic <- function(x, params_0, eps = 0.001)
{
  params <- params_0
  n <- length(x)
  
  sumlogx <- sum(log(x))
  
  diff <- Inf
  alpha <- params_0[1]
  beta <- params_0[2]
  
  while (diff > eps)
  {
    params.old <- params
    
    # Intermediary variables
    w1 <- sum((x/alpha)^beta * (1 + (x/alpha)^beta)^-1)
    w2 <- sum(log(x) * (x/alpha)^beta * (1 + (x/alpha)^beta)^-1)
    w3 <- sum(log(x)^2 * (x/alpha)^beta * (1 + (x/alpha)^beta)^-1)
    
    # Gradient vector
    logx_over_alpha <- log(x/alpha)
    if (any(is.nan(logx_over_alpha))) {
      s <- c(0, 0)
    } else {
      s <- c(n/beta - sum(log(1 + (x/alpha)^beta)), 
             -n*alpha/beta^2 +
               sum(logx_over_alpha * (x/alpha)^beta * (1 + (x/alpha)^beta)^-2))
    }
    
    # Jacobian matrix
    if (any(is.nan(logx_over_alpha))) {
      J <- matrix(0, nrow=2, ncol=2)
    } else {
      J <- matrix(0, nrow=2, ncol=2)
      J[1, 1] <- n/alpha^2 + beta * sum((x/alpha)^(2*beta)/(1 + (x/alpha)^beta)^2)
      J[1, 2] <- -sum((x/alpha)^beta * logx_over_alpha/(1 + (x/alpha)^beta))
      J[2, 1] <- -beta * sum((x/alpha)^beta * logx_over_alpha/(1 + (x/alpha)^beta))
      J[2, 2] <- n/beta^2 - sum((log(1 + (x/alpha)^beta))^2)
    }
    
    # Update parametres using jacobian matrix
    params <- params - solve(J) %*% s
    
    alpha <- params[1]
    beta <- params[2]
    
    diff <- abs(sum(params - params.old))
  }
  
  list(params = params, J = J)
}


gpigs <- read.table("surv.gpigs.txt", header = T, sep = ";")
gpigs.m43 <- gpigs.noncensored$lifetime[gpigs.noncensored$regime == "M_4.3"]

NR_fit_LogLogistic(gpigs.m43, c(150, 3))
# pour comparaison le résultat de fitdist:  alpha = 152.389696   beta = 3.012415 
# on est bien loin

# next try avec la matrice hessienne
loglogistic_hessian <- function(x, n, alpha, beta) {
  hessian <- matrix(0, nrow = 2, ncol = 2)
  for (i in 1:length(x)) {
    xi <- x[i]
    hessian[1, 1] <- hessian[1, 1] + (beta^2 / xi^2) * ((xi / alpha)^(beta - 1)) / ((1 + (xi / alpha)^beta)^2)
    hessian[1, 2] <- hessian[1, 2] + (2 * xi^beta * log(xi / alpha)) / ((1 + (xi / alpha)^beta)^2)
    hessian[2, 1] <- hessian[2, 1] + (2 * xi^beta * log(xi / alpha)) / ((1 + (xi / alpha)^beta)^2)
    hessian[2, 2] <- hessian[2, 2] - (xi^beta * log(xi / alpha)) / ((1 + (xi / alpha)^beta)^2)
  }
  hessian[1, 1] <- (n * beta / alpha^2) - (2 * hessian[1, 1])
  hessian[2, 2] <- -(n / beta^2) - (2 * hessian[2, 2])
  hessian[2, 1] <- -(n / alpha) + (2 * hessian[2, 1])
  hessian[1, 2] <- hessian[2, 1]
  hessian
}


NR_fit_LogLogistic_hessian <- function(x, params_0, eps = 0.001)
{
  params <- params_0
  n <- length(x)
  
  sumlogx <- sum(log(x))
  
  diff <- Inf
  alpha <- params_0[1]
  beta <- params_0[2]
  
  while (diff > eps)
  {
    params.old <- params
    
    # Intermediary variables
    w1 <- sum((x/alpha)^beta * (1 + (x/alpha)^beta)^-1)
    w2 <- sum(log(x) * (x/alpha)^beta * (1 + (x/alpha)^beta)^-1)
    w3 <- sum(log(x)^2 * (x/alpha)^beta * (1 + (x/alpha)^beta)^-1)
    
    # Gradient vector
    logx_over_alpha <- log(x/alpha)
    if (any(is.nan(logx_over_alpha))) {
      s <- c(0, 0)
    } else {
      s <- c(n/beta - sum(log(1 + (x/alpha)^beta)), 
             -n*alpha/beta^2 +
               sum(logx_over_alpha * (x/alpha)^beta * (1 + (x/alpha)^beta)^-2))
    }
    
    # Hessian matrix
    H <- loglogistic_hessian(x, n, alpha, beta)
    
    # Update parametres using jacobian matrix
    params <- params - solve(H) %*% s
    
    alpha <- params[1]
    beta <- params[2]
    
    diff <- abs(sum(params - params.old))
  }
  
  list(params = params, H = H)
}

NR_fit_LogLogistic_hessian(gpigs.m43, c(150, 3))
# IT WOOOOORKS WOOHOOOOO
