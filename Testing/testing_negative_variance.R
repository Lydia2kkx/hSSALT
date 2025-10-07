library(hSSALT)

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
seed <- 1:1000



############Continuous Type-I case
results.cont.asymptotic <- list()
results.cont.percentile <- list()
results.cont.bca <- list()

MLE_NotNA.cont <- rep(0, length(seed))
MLE_converge.cont <- rep(0, length(seed))
MLE_heter.cont <- rep(0, length(seed))

for (i in seed) {
  set.seed(i)
  sample <- rhSSALT(n,1,tau=tau,theta1=theta1,theta21 = theta21,theta22 = theta22,p=p, monitoring="continuous")
  resMLE <- MLEhSSALT(sample$Censored_dat,n,1,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")
  
  # Test the Sample
  if (!is.na(resMLE$info)) {
    MLE_NotNA.cont[i] <- TRUE
    MLE_heter.cont[i] <- ifelse(resMLE$info == "heterogeneous", TRUE, FALSE)
    MLE_converge.cont[i] <- ifelse(resMLE$message == "convergent", TRUE, FALSE)
  } else {
    MLE_NotNA.cont[i] <- FALSE
    MLE_converge.cont[i] <- FALSE
    MLE_heter.cont[i] <- FALSE
  }
  
  #if (is.na(resMLE$info)  || resMLE$info != "heterogeneous" || resMLE$message != "convergent") {
  if (is.na(resMLE$info)) {
    results.cont.asymptotic[[i]] <- list(theta1=c(0,0), theta21=c(0,0), theta22=c(0,0),p=c(0,0), B=0, j=0, alpha=0, type=0)
    results.cont.percentile[[i]] <- list(theta1=c(0,0), theta21=c(0,0), theta22=c(0,0),p=c(0,0), B=0, j=0, alpha=0, type=0)
    results.cont.bca[[i]] <- list(theta1=c(0,0), theta21=c(0,0), theta22=c(0,0),p=c(0,0), B=0, j=0, alpha=0, type=0)
  } else {
    
    ### Results
    #Asymptotic
    
    set.seed(i)
    sample <- rhSSALT(n,1,tau=tau,theta1=theta1,theta21 = theta21,theta22 = theta22,p=p, monitoring="continuous")
    resMLE <- MLEhSSALT(sample$Censored_dat,n,1,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")
    
    results.cont.asymptotic[[i]] <- tryCatch(
      {
        CIhSSALT(sample$Censored_dat,n,1,tau=tau,MLEhSSALT_Obj=resMLE,language = "CPP",monitoring = "continuous", B=B,grid=F)
      },
      error = function(e) {
        list(theta1=c(NA,NA), theta21=c(NA,NA), theta22=c(NA,NA),p=c(NA,NA), B=NA, j=NA, alpha=NA, type=NA)
      }
    )
  }
  cat("i: ", i, "\n")
}

### End Results

#Asymptotic
df.cont.asymptotic.temp <- do.call(rbind, lapply(results.cont.asymptotic, function(x) {
  data.frame(
    theta1_l = x$theta1[1],
    theta1_u = x$theta1[2],
    theta21_l = x$theta21[1],
    theta21_u = x$theta21[2],
    theta22_l = x$theta22[1],
    theta22_u = x$theta22[2],
    p_l = x$p[1],
    p_u = x$p[2],
    B=x$B,
    j=x$j,
    alpha=x$alpha,
    type=x$type
  )
}))
df.cont.asymptotic.cpp <- cbind(seed, df.cont.asymptotic.temp, MLE_NotNA.cont, MLE_converge.cont, MLE_heter.cont)
df.cont.asymptotic.cpp

# ############Interval case
# results.int.asymptotic <- list()
# results.int.percentile <- list()
# results.int.bca <- list()
# 
# MLE_NotNA.int <- rep(0, length(seed))
# MLE_converge.int <- rep(0, length(seed))
# MLE_heter.int <- rep(0, length(seed))
# 
# for (i in seed) {
#   set.seed(i)
#   sample <- rhSSALT(n,1,tau=tau,theta1=theta1,theta21 = theta21,theta22 = theta22,p=p, monitoring="interval",delta=delta)
#   resMLE <- MLEhSSALT(sample$Censored_dat,n,1,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="interval", delta=delta)
#   
#   # Test the Sample
#   if (!is.na(resMLE$info)) {
#     MLE_NotNA.int[i] <- TRUE
#     MLE_heter.int[i] <- ifelse(resMLE$info == "heterogeneous", TRUE, FALSE)
#     MLE_converge.int[i] <- ifelse(resMLE$message == "convergent", TRUE, FALSE)
#   } else {
#     MLE_NotNA.int[i] <- FALSE
#     MLE_converge.int[i] <- FALSE
#     MLE_heter.int[i] <- FALSE
#     
#   }
#   
#   if (is.na(resMLE$info)  || resMLE$info != "heterogeneous" || resMLE$message != "convergent") {
#     results.int.asymptotic[[i]] <- list(theta1=c(0,0), theta21=c(0,0), theta22=c(0,0),p=c(0,0), B=0, j=0, alpha=0, type=0)
#     results.int.percentile[[i]] <- list(theta1=c(0,0), theta21=c(0,0), theta22=c(0,0),p=c(0,0), B=0, j=0, alpha=0, type=0)
#     results.int.bca[[i]] <- list(theta1=c(0,0), theta21=c(0,0), theta22=c(0,0),p=c(0,0), B=0, j=0, alpha=0, type=0)
#   } else {
#     
#     ### Results
#     #Asymptotic
#     
#     set.seed(i)
#     sample <- rhSSALT(n,1,tau=tau,theta1=theta1,theta21 = theta21,theta22 = theta22,p=p, monitoring="interval",delta=delta)
#     resMLE <- MLEhSSALT(sample$Censored_dat,n,1,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="interval", delta=delta)
#     
#     results.int.asymptotic[[i]] <- tryCatch(
#       {
#         CIhSSALT(sample$Censored_dat,n,1,tau=tau,MLEhSSALT_Obj=resMLE,language = "CPP",monitoring = "interval",delta=delta, B=B,grid=F)
#       },
#       error = function(e) {
#         list(theta1=c(NA,NA), theta21=c(NA,NA), theta22=c(NA,NA),p=c(NA,NA), B=NA, j=NA, alpha=NA, type=NA)
#       }
#     )
#   }
#   cat("i: ", i, "\n")
# }
# 
# ### End Results
# 
# #Asymptotic
# df.int.asymptotic.temp <- do.call(rbind, lapply(results.int.asymptotic, function(x) {
#   data.frame(
#     theta1_l = x$theta1[1],
#     theta1_u = x$theta1[2],
#     theta21_l = x$theta21[1],
#     theta21_u = x$theta21[2],
#     theta22_l = x$theta22[1],
#     theta22_u = x$theta22[2],
#     p_l = x$p[1],
#     p_u = x$p[2],
#     B=x$B,
#     j=x$j,
#     alpha=x$alpha,
#     type=x$type
#   )
# }))
# df.int.asymptotic.cpp <- cbind(seed, df.int.asymptotic.temp, MLE_NotNA.int, MLE_converge.int, MLE_heter.int)
# df.int.asymptotic.cpp