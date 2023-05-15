
#####################################
###          MONTECARLO           ###
#####################################
source("fitdistrmachin.r")
source("newton-raphson.r")
set.seed(42)
#####################################
# some values needed for later
n_values <- c(100, 500, 1000, 5000, 10000) # sample sizes
I <- 1000 # no. of iterations
alpha0 <- 452.6
beta0 <- 2.2
#####################################
# censoring function, for a given quantile sets the max value and censors everything above
censor <- function(x, p = 0.9) {
  # for example the top 10% is everything above the 0.9 quantile
  top_ten <- quantile(x, p)
  res = x
  # censor that shit
  res[x > top_ten] = top_ten
  return(res)
}
#####################################
# the montecarlo iteration for a given sample size
mc <- function(n, I, alpha0, beta0) {
  # initialise the matrices that will track the estimates throughout the iterations
  est <- matrix(nrow = 0, ncol = 2)
  # iterate
  for (i in 1:I) {
    print(c(i, I))
    # create random data set that follows a log-logistic distribution with the given parametres
    dat <- rloglogis(n, alpha0, beta0)
    # censor the top 10%
    dat_censored <- censor(dat)
    est_temp <- NR_fit_LogLogistic_hessian(dat_censored, c(alpha0, beta0), 0.0001)$params
    est <- rbind(est, c(est_temp[1,1], est_temp[2,1]))
  }
  return(c(mean((est[,1] - alpha0)^2), mean((est[,2] - beta0)^2)))
}
# test_mc <- mc(100, I, alpha0, beta0)
#####################################
# res <- matrix(nrow = 2, ncol = 0)
# for (i in 1:length(n_values)) {
#   res <- cbind(res, mc(n_values[i], I, alpha0, beta0))
# }
# colnames(res) <- n_values[1:4]
# write.table(res, "mc_res_with_censoring_1k_iterations_452_2.csv")
#####################################
res <- as.matrix(read.table("mc_res_with_censoring_1k_iterations_452_2.csv",header=T))

barplot(log(res),
        col = c("blue", "orange"),
        beside = TRUE,
        main = "Mean square error of the NR-fit parameter estimates \n based on sample size for the \n Log-Logistic distribution \n (Logarithmic scale)", xlab = "Sample size", ylab = "log(MSE)")
legend("topright",
       legend = c("Alpha", "Beta"),
       fill = c("blue", "orange"))
