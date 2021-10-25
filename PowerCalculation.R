simFit <- function(n){
  SP <- CNVDATA$Begin
  EP <- CNVDATA$End
  X1 <- EP-SP
  Type <- CNVDATA$Type
  X2 <- Type
  Dosage <- CNVDATA$Dosage
  X3 <- abs(Dodage-2)
  xb <- Beta0 + Beta*X1*X2*X3
  p <- 1/(1 + exp(-xb))
  y <- rbinom(n = n, size = 1, prob = p)
  mod <- glm(y ~ X1 + X2 + X3 + X1:X2:X3, family = binomial)
  s.out <- summary(mod)
  s.out$coefficients["X1:X2:X3
                     ","Pr(>|z|)"] < 0.05
}

powerEst <- function(N, n){
  r.out <- replicate(n = N, simFit(n = n))
  mean(r.out)
}

powerEst(N = 500, n = 200)

ss <- seq(300,1000,100) # various sample sizes
p.out <- sapply(ss, function(x)powerEst(N = 500, n = x))

plot(ss, p.out, type = "l", 
     xlab = "sample size", ylab = "power")
abline(h = 0.8, lty = 2)
