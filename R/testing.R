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


############Continuous Type-I case
results <- list()
m <- 1

for (i in 1:100) {
  set.seed(i)
  n <- 40
  #sample <- rhSSALT(n,1,tau=c(8,16),theta1=33.12,theta21 = 2,theta22 = 15,p=0.4, monitoring="interval",delta=0.5)
  sample <- rhSSALT(n,1,tau=c(8,16),theta1=33.12,theta21 = 2,theta22 = 15,p=0.4)
  #sample <- rhSSALT(n,2,r=30,tau=8,theta1=33,theta21 = 2,theta22 = 7.4,p=0.4)
  #print(sample)
  
  #resMLE <- MLEhSSALT(sample$Censored_dat,n,1,tau=c(8,16),theta21 = 2,theta22 = 15,p=0.4,language = "R", monitoring="interval", delta=0.5)
  resMLE <- MLEhSSALT(sample$Censored_dat,n,1,tau=c(8,16),theta21 = 2,theta22 = 15,p=0.4,language = "R", monitoring="continuous")
  #resMLE <- MLEhSSALT(sample$Censored_dat,n,2,r=30,tau=8,theta21 = 2,theta22 = 7,p=0.4,language = "R", monitoring="continuous")
  #print(resMLE)
  
  if (!is.na(resMLE$info)) {
    if (resMLE$info == "heterogeneous" && resMLE$message == "convergent") {
      cat("i: ", i, "\n")
      #res <- CIhSSALT(sample$Censored_dat,n,1,tau=c(8,16),MLEhSSALT_Obj=resMLE,language = "R",monitoring = "interval",delta=0.5, B=500,CImethod="bca",grid=FALSE)
      res <- CIhSSALT(sample$Censored_dat,n,1,tau=c(8,16),MLEhSSALT_Obj=resMLE,language = "R", CImethod="bca", B=500,grid=FALSE)
      results[[m]] <- res
      m <- m+1
    } 
  }
  #res <- CIhSSALT(sample$Censored_dat,n,1,tau=c(8,16),MLEhSSALT_Obj=resMLE,language = "R",monitoring = "interval",delta=1, B=100,CImethod="bca")
  #res <- CIhSSALT(sample$Censored_dat,n,1,tau=c(8,16),MLEhSSALT_Obj=resMLE,language = "R", CImethod="bca", B=500)
  #res <- CIhSSALT(sample$Censored_dat,n=n,MLEhSSALT_Obj=resMLE,r=30,censoring=2,tau=8,language = "R", CImethod="bca", B=100)
  #print(res)
  #i <- i+1
  #cat("i: ", i, "\n")
}


############Interval case
results <- list()
m <- 1

for (i in 1:100) {
  set.seed(i)
  n <- 40
  sample <- rhSSALT(n,1,tau=c(8,16),theta1=33.12,theta21 = 2,theta22 = 15,p=0.4, monitoring="interval",delta=0.5)
  #sample <- rhSSALT(n,1,tau=c(8,16),theta1=33.12,theta21 = 2,theta22 = 15,p=0.4)
  #sample <- rhSSALT(n,2,r=30,tau=8,theta1=33,theta21 = 2,theta22 = 7.4,p=0.4)
  #print(sample)
  
  resMLE <- MLEhSSALT(sample$Censored_dat,n,1,tau=c(8,16),theta21 = 2,theta22 = 15,p=0.4,language = "R", monitoring="interval", delta=0.5)
  #resMLE <- MLEhSSALT(sample$Censored_dat,n,1,tau=c(8,16),theta21 = 2,theta22 = 15,p=0.4,language = "R", monitoring="continuous")
  #resMLE <- MLEhSSALT(sample$Censored_dat,n,2,r=30,tau=8,theta21 = 2,theta22 = 7,p=0.4,language = "R", monitoring="continuous")
  #print(resMLE)
  
  if (!is.na(resMLE$info)) {
    if (resMLE$info == "heterogeneous" && resMLE$message == "convergent") {
      cat("i: ", i, "\n")
      res <- CIhSSALT(sample$Censored_dat,n,1,tau=c(8,16),MLEhSSALT_Obj=resMLE,language = "R",monitoring = "interval",delta=0.5, B=500,CImethod="bca",grid=FALSE)
      #res <- CIhSSALT(sample$Censored_dat,n,1,tau=c(8,16),MLEhSSALT_Obj=resMLE,language = "R", CImethod="bca", B=500,grid=FALSE)
      results[[m]] <- res
      m <- m+1
    } 
  }
  #res <- CIhSSALT(sample$Censored_dat,n,1,tau=c(8,16),MLEhSSALT_Obj=resMLE,language = "R",monitoring = "interval",delta=1, B=100,CImethod="bca")
  #res <- CIhSSALT(sample$Censored_dat,n,1,tau=c(8,16),MLEhSSALT_Obj=resMLE,language = "R", CImethod="bca", B=500)
  #res <- CIhSSALT(sample$Censored_dat,n=n,MLEhSSALT_Obj=resMLE,r=30,censoring=2,tau=8,language = "R", CImethod="bca", B=100)
  #print(res)
  #i <- i+1
  #cat("i: ", i, "\n")
}

