gpigs <- read.table("surv.gpigs.txt", header = T, sep = ";")
gpigs.noncensored <- gpigs[gpigs$censored == 0,]

gpigs.mc <- gpigs.noncensored$lifetime[gpigs.noncensored$regime == "M_C"]
gpigs.mc.with_censored <- gpigs$lifetime[gpigs$regime == "M_C"]
gpigs.m43 <- gpigs.noncensored$lifetime[gpigs.noncensored$regime == "M_4.3"]

summary(gpigs.mc)
summary(gpigs.mc.with_censored) # on garde celui-là pour la summary
summary(gpigs.m43)

hist(gpigs.noncensored$lifetime, breaks = 26)

dloglogis <- function(x, alpha, beta) {
  res <- (beta/alpha) * (x/alpha)^(beta-1)*(1+(x/alpha)^beta)^-2
  return(res)
}

ploglogis <- function(q, alpha, beta) {
  res <- 1 / (1 + (q/alpha)^beta)
  return(res)
}

qloglogis <- function(p, alpha, beta) {
  res <- alpha * ((1/p) - 1)^(1/beta)
  return(res)
}


library(fitdistrplus)
# MC regime (sans données censorées pour le dist)
fw <- fitdist(gpigs.mc, "loglogis", method = "mle", start=list(alpha=10, beta=5))
fw$estimate
denscomp(list(fw), main = "MC regime")

# M4.3 regime
fw <- fitdist(gpigs.noncensored$lifetime[gpigs.noncensored$regime == "M_4.3"], "loglogis", method = "mle", start=list(alpha=10, beta=5))
fw$estimate
denscomp(list(fw), main = "M4.3 regime")
