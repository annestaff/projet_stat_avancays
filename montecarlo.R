
#####################################
###          MONTECARLO           ###
#####################################
source("fitdistrmachin.r")
source("newton-raphson.r")
set.seed(42)
#####################################
n_values <- c(100, 500, 1000, 5000, 10000)
I <- 1000
alpha0 <- 300
beta0 <- 5
#####################################
mc <- function(n, I, alpha0, beta0) {
  est <- matrix(nrow = 0, ncol = 2)
  mse <- matrix(nrow = 0, ncol = 2)
  for (i in 1:I) {
    dat <- rloglogis(n, alpha0, beta0)
    est_temp <- NR_fit_LogLogistic_hessian(dat, c(alpha0, beta0))$params
    est <- rbind(est, c(est_temp[1,1], est_temp[2,1]))
    mse_alpha <- mean((est[i,1] - alpha0)^2)
    mse_beta <- mean((est[i,2] - beta0)^2)
    mse <- rbind(mse, c(mse_alpha, mse_beta))
  }
  return(colMeans(mse))
}
test_mc <- mc(n, I, alpha_0, beta_0)
#####################################
res <- matrix(nrow = 2, ncol = 0)
for (i in 1:length(n_values)) {
  res <- cbind(res, mc(n_values[i], I, alpha0, beta0))
}
colnames(res) <- n_values
write.table(res, "mc_res.csv")
#####################################
# res <- read.table("mc_res.csv")

barplot(log(res),
        col = c("blue", "orange"),
        beside = TRUE,
        main = "Mean square error of the NR-fit parameter estimates \n based on sample size for the \n Log-Logistic distribution \n (Logarithmic scale)", xlab = "Sample size", ylab = "log(MSE")
legend("topright",
       legend = c("Alpha", "Beta"),
       fill = c("blue", "orange"))