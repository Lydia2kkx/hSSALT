#boot
library(boot)

x <- rnorm(100, mean = 10, sd = 2)

boot_mean <- function(data, ind) {
  mean(data[ind])
}

b <- boot(data = x, statistic = boot_mean, R = 1000)
boot.ci(b, conf = 0.95, type = c("perc", "bca"))


#car
library(car)

x1 <- runif(100)
x2 <- runif(100)
y  <- 1 + 2*x1 - 3*x2 + rnorm(100, sd = 0.3)
regmod <- lm(y ~ x1 + x2)

boot_mod <- Boot(regmod, R = 1000, method = "case")
confint(boot_mod, type = "bca")

#simpleboot
library(simpleboot)

x <- runif(100)

boot_res <- one.boot(x, mean, R = 1000)
boot.ci(boot_res, type = "perc")


#hSSALT
library(hSSALT)

n <- 40
tau <- c(8,16)
theta1 <- 33.12
theta21 <- 2
theta22 <- 15
p <- 0.4

set.seed(3)
sample <- rhSSALT(n,1,tau=tau,theta1=theta1,theta21 = theta21,theta22 = theta22,p=p, monitoring="continuous")
resMLE_cpp <- MLEhSSALT(sample$Censored_dat,n,1,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")
CI_cpp <- CIhSSALT(sample$Censored_dat,n,1,tau=tau,MLEhSSALT_Obj=resMLE_cpp,language = "CPP",monitoring="continuous", B=500, CImethod="percentile")
CI_cpp