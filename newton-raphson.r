### NEWTON RAPHSON ALGORITHM ###
# Goal: find the root (x axis intersection) of a function, works for real and complex functions
# 1. pick a random point x of the function ("first guess")
# 2. find the tangent to the curve at that point x
# 3. find the intersection of the tangent with the x axis and make it your new x
# 4. repeat until you're close enough (delta-y below a certain threshold which defines that you're on a point)

score_alpha <- function(x, alpha, beta) {
  n <- length(x)
  score <- -((n * beta) / alpha) + 2 * sum((beta / x) * ((x / alpha)^(beta - 1)) / (1 + (x / alpha)^beta))
  return(score)
}

score_beta <- function(x, alpha, beta) {
  n <- length(x)
  score <- (n / beta) - (n * log(alpha)) + sum(log(x)) + 2 * sum(((x / alpha)^beta) * log(x / alpha) / (1 + (x / alpha)^beta))
  return(score)
}

loglogistic_hessian <- function(x, n, alpha, beta) {
  hessian <- matrix(0, nrow = 2, ncol = 2)
  xi <- x
  hessian[1, 1] <- (n * beta / alpha^2) - (2 * sum((beta^2 / xi^2) * ((xi / alpha)^(beta - 1)) / ((1 + (xi / alpha)^beta)^2)))
  hessian[2, 2] <- -(n / beta^2) - (2 * sum(- (xi^beta * log(xi / alpha)) / ((1 + (xi / alpha)^beta)^2)))
  hessian[2, 1] <- -(n / alpha) + (2 * sum((2 * xi^beta * log(xi / alpha)) / ((1 + (xi / alpha)^beta)^2)))
  hessian[1, 2] <- sum((2 * xi^beta * log(xi / alpha)) / ((1 + (xi / alpha)^beta)^2))
  return(hessian)
}


NR_fit_LogLogistic_hessian <- function(x, params_0, eps = 0.00001)
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
    # Gradient vector
    s <- c(score_alpha(x, alpha, beta), score_beta(x, alpha, beta))
    
    # Calculate Hessian matrix
    H <- loglogistic_hessian(x, n, alpha, beta)
    
    # Update parametres using hessian matrix
    params <- params - solve(H) %*% s
    alpha <- params[1]
    beta <- params[2]
    
    diff <- abs(sum(params - params.old))
  }
  
  list(params = params, H = H)
}

# NR_fit_LogLogistic_hessian(gpigs.m43, c(150, 3))
# pour comparaison le résultat de fitdist:  alpha = 152.389696   beta = 3.012415 
