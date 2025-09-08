source("rhSSALT.R")
source("MLEhSSALT.R")
source("MLE_Exp.R")
source("EM_algorithm_censored.R")
source("EM_algorithm_interval.R")
source("MLE_Geo.R")
source("bs_bca_continuous_help_functions.R")
source("CIbca_hSSALT.R")
source("CIbs_hSSALT.R")
source("CIhSSALT.R")


############Variables
n <- 40
tau <- c(8,16)
theta1 <- 33.12
theta21 <- 2
theta22 <- 15
p <- 0.4
B <- 500
r <- 30
delta <- 0.5

############Continuous Type-I case
results <- list()
m <- 1

for (i in 1:100) {
  set.seed(i)
  sample <- rhSSALT(n,1,tau=tau,theta1=theta1,theta21 = theta21,theta22 = theta22,p=p)
  resMLE <- MLEhSSALT(sample$Censored_dat,n,1,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "R", monitoring="continuous")

  if (!is.na(resMLE$info)) {
    if (resMLE$info == "heterogeneous" && resMLE$message == "convergent") {
      cat("i: ", i, "\n")
      res <- CIhSSALT(sample$Censored_dat,n,1,tau=tau,MLEhSSALT_Obj=resMLE,language = "R", CImethod="bca", B=B,grid=T)
      results[[m]] <- res
      m <- m+1
    }
  }
}

print(results)

############Continuous Type-II case
results <- list()
m <- 1

for (i in 1:100) {
  set.seed(i)
  sample <- rhSSALT(n,2,r=r,tau=tau[1],theta1=theta1,theta21 = theta21,theta22 = theta22,p=p)
  resMLE <- MLEhSSALT(sample$Censored_dat,n,2,r=r,tau=tau[1],theta21 = theta21,theta22 = theta22,p=p,language = "R", monitoring="continuous")
  
  if (!is.na(resMLE$info)) {
    if (resMLE$info == "heterogeneous" && resMLE$message == "convergent") {
      cat("i: ", i, "\n")
      res <- CIhSSALT(sample$Censored_dat,n,1,tau=c(8,16),MLEhSSALT_Obj=resMLE,language = "R", CImethod="bca", B=500,grid=T) #Modify for Type-II
      results[[m]] <- res
      m <- m+1
    }
  }
}

print(results)


############Interval case
results <- list()
m <- 1

for (i in 1:100) {
  set.seed(i)
  sample <- rhSSALT(n,1,tau=tau,theta1=theta1,theta21 = theta21,theta22 = theta22,p=p, monitoring="interval",delta=delta)
  resMLE <- MLEhSSALT(sample$Censored_dat,n,1,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "R", monitoring="interval", delta=delta)

  if (!is.na(resMLE$info)) {
    if (resMLE$info == "heterogeneous" && resMLE$message == "convergent") {
      cat("i: ", i, "\n")
      res <- CIhSSALT(sample$Censored_dat,n,1,tau=tau,MLEhSSALT_Obj=resMLE,language = "R",monitoring = "interval",delta=delta, B=B,CImethod="bca",grid=F)
      results[[m]] <- res
      m <- m+1
    }
  }
}

print(results)
