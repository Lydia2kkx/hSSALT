source("../R/rhSSALT.R")
source("../R/MLEhSSALT.R")
source("../R/MLE_Exp.R")
source("../R/EM_algorithm_censored.R")
source("../R/EM_algorithm_interval.R")
source("../R/MLE_Geo.R")
source("../R/bs_bca_continuous_help_functions.R")
source("../R/CIbca_hSSALT.R")
source("../R/CIbs_hSSALT.R")
source("../R/CIsay_hSSALT.R")
source("../R/CIhSSALT.R")
source("../R/print.CIhSSALT.R")
source("../R/print.hSSALTMLE.R")

library(Rcpp)

sourceCpp("../src/EM_arma_int_cont.cpp")


############Variables from PDF
n <- 35
tau <- 8
theta1 <- exp(3.5)
theta21 <- exp(-0.2)
theta22 <- exp(2)
p <- 0.4
B <- 1000
r <- 30
#delta <- 0.5

data <- c(1.46,2.22,3.92,4.24,5.47,5.6,6.12,6.56,8.01,8.1,8.22,8.3,8.59,8.77,8.8,8.8,8.84,8.9,8.97,8.98,9.43,9.62,
          9.87,11.14,11.85,12.14,13.49,14.19,14.33,15.28,16.58,17.8,21.09,26.34,28.66)


#########Type I - Continuous
#### Given results Dataframe
th1_l <- c(8.43,18.91,18.95,4.93,17.56,17.5,0,15.21,14.43)
th1_u <- c(45.01,65.62,65.68,48.51,67.32,67.13,55.36,91,90.3)
th21_l <- c(0.23,0.29,0.3,0.12,0.21,0.22,0,0.06,0.06)
th21_u <- c(1.43,1.8,1.86,1.54,2.03,2.06,1.77,2.92,2.92)
th22_l <- c(3.8,4.33,4.38,3.13,3.9,3.9,1.82,3.23,2.75)
th22_u <- c(10.82,14.87,15.34,11.5,17.56,17.55,12.81,25.14,23.54)
p_l <- c(0.18,0.22,0.22,0.14,0.19,0.19,0.06,0.12,0.11)
p_u <- c(0.62,0.67,0.68,0.67,0.71,0.72,0.75,0.77,0.8)

alpha_c <- c(rep(0.1,3),rep(0.05,3),rep(0.01,3))
type_c <- rep(c("Asymptotic", "Percentile", "BCa"), 3)

df.i <- data.frame(
  th1_l = th1_l,
  th1_u = th1_u,
  th21_l = th21_l,
  th21_u = th21_u,
  th22_l = th22_l,
  th22_u = th22_u,
  p_l = p_l,
  p_u = p_u,
  alpha_c = alpha_c,
  type_c = type_c
)

# resMLE.3 <- MLEhSSALT(data,n,1,tau=c(8,16),theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")
# resMLE.4 <- MLEhSSALT(data,n,1,tau=c(8,20),theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")

#alpha <- 0.90
resMLE.5 <- MLEhSSALT(data,n,1,tau=c(8,24),theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")

CI_i.11 <- CIhSSALT(data, n, 1, tau = c(8,24), MLEhSSALT_Obj = resMLE.5, language = "CPP",monitoring = "continuous", B = B, alpha = 0.1)
row1 <- round(as.numeric(unlist(CI_i.11)[-c(9,10,12)]),2)

CI_i.12 <- CIhSSALT(data, n, 1, tau = c(8,24), MLEhSSALT_Obj = resMLE.5, language = "CPP",monitoring = "continuous", B = B, CImethod="percentile", alpha = 0.1)
row2 <- round(as.numeric(unlist(CI_i.12)[-c(9,10,12)]),2)

CI_i.13 <- CIhSSALT(data, n, 1, tau = c(8,24), MLEhSSALT_Obj = resMLE.5, language = "CPP",monitoring = "continuous", B = B, CImethod="bca", alpha = 0.1)
row3 <- round(as.numeric(unlist(CI_i.13)[-c(9,10,12)]),2)

#alpha <- 0.95
resMLE.5 <- MLEhSSALT(data,n,1,tau=c(8,24),theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")

CI_i.21 <- CIhSSALT(data, n, 1, tau = c(8,24), MLEhSSALT_Obj = resMLE.5, language = "CPP",monitoring = "continuous", B = B)
row4 <- round(as.numeric(unlist(CI_i.21)[-c(9,10,12)]),2)

CI_i.22 <- CIhSSALT(data, n, 1, tau = c(8,24), MLEhSSALT_Obj = resMLE.5, language = "CPP",monitoring = "continuous", B = B, CImethod="percentile")
row5 <- round(as.numeric(unlist(CI_i.22)[-c(9,10,12)]),2)

CI_i.23 <- CIhSSALT(data, n, 1, tau = c(8,24), MLEhSSALT_Obj = resMLE.5, language = "CPP",monitoring = "continuous", B = B, CImethod="bca")
row6 <- round(as.numeric(unlist(CI_i.23)[-c(9,10,12)]),2)

#alpha <- 0.99
resMLE.5 <- MLEhSSALT(data,n,1,tau=c(8,24),theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")

CI_i.31 <- CIhSSALT(data, n, 1, tau = c(8,24), MLEhSSALT_Obj = resMLE.5, language = "CPP",monitoring = "continuous", B = B, alpha = 0.01)
row7 <- round(as.numeric(unlist(CI_i.31)[-c(9,10,12)]),2)

CI_i.32 <- CIhSSALT(data, n, 1, tau = c(8,24), MLEhSSALT_Obj = resMLE.5, language = "CPP",monitoring = "continuous", B = B, CImethod="percentile", alpha = 0.01)
row8 <- round(as.numeric(unlist(CI_i.32)[-c(9,10,12)]),2)

CI_i.33 <- CIhSSALT(data, n, 1, tau = c(8,24), MLEhSSALT_Obj = resMLE.5, language = "CPP",monitoring = "continuous", B = B, CImethod="bca", alpha = 0.01)
row9 <- round(as.numeric(unlist(CI_i.33)[-c(9,10,12)]),2)

df.i.funcs <- rbind(row1, row2, row3, row4, row5, row6, row7, row8, row9)
df.i.funcs <- as.data.frame(df.i.funcs, stringsAsFactors = FALSE)
colnames(df.i.funcs) <- colnames(df.i)[1:9]
rownames(df.i.funcs) <- NULL
df.i.funcs$type_c <- type_c

#Make Difference Dataframe
df.i.diffs <- df.i
df.i.diffs[1:8] <- round(abs(df.i[1:8] - df.i.funcs[1:8]), 2)
df.i.diffs



#############Type II - r=n
#### Given results Dataframe
th1_l <- c(8.41,18.86,18.77,4.90,17.48,17.37,0,14.9,12.94)
th1_u <- c(44.98,66.22,64.65,48.49,68.28,66.83,55.34,92.63,90.53)
th21_l <- c(0.16,0.19,0.24,0.04,0.1,0.14,0,0.01,0.01)
th21_u <- c(1.43,1.84,2.02,1.55,2.26,2.53,1.79,3.65,3.78)
th22_l <- c(3.34,4.1,4.48,2.71,3.72,4.11,1.49,2.81,3.61)
th22_u <- c(9.89,11.52,14.3,10.51,14.29,17.23,11.74,18.71,24.76)
p_l <-  c(0.11,0.14,0.14,0.06,0.09,0.11,0,0.06,0.06)
p_u <- c(0.65,0.69,0.73,0.7,0.8,0.82,0.8,0.91,0.91)

alpha_c <- c(rep(0.1,3),rep(0.05,3),rep(0.01,3))
type_c <- rep(c("Asymptotic", "Percentile", "BCa"), 3)

df.ii <- data.frame(
  th1_l = th1_l,
  th1_u = th1_u,
  th21_l = th21_l,
  th21_u = th21_u,
  th22_l = th22_l,
  th22_u = th22_u,
  p_l = p_l,
  p_u = p_u,
  alpha_c = alpha_c,
  type_c = type_c
)

#alpha <- 0.90
resMLE.1 <- MLEhSSALT(data,n,2,r=n,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")

CI_i.11 <- CIhSSALT(data, n, 2, r=n, tau = tau, MLEhSSALT_Obj = resMLE.1, language = "CPP",monitoring = "continuous", B = B, alpha = 0.1)
row1 <- round(as.numeric(unlist(CI_i.11)[-c(9,10,12)]),2)

CI_i.12 <- CIhSSALT(data, n, 2, r=n, tau = tau, MLEhSSALT_Obj = resMLE.1, language = "CPP",monitoring = "continuous", B = B, CImethod="percentile", alpha = 0.1)
row2 <- round(as.numeric(unlist(CI_i.12)[-c(9,10,12)]),2)

CI_i.13 <- CIhSSALT(data, n, 2, r=n, tau = tau, MLEhSSALT_Obj = resMLE.1, language = "CPP",monitoring = "continuous", B = B, CImethod="bca", alpha = 0.1)
row3 <- round(as.numeric(unlist(CI_i.13)[-c(9,10,12)]),2)

#alpha <- 0.95
resMLE.1 <- MLEhSSALT(data,n,2,r=n,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")

CI_i.21 <- CIhSSALT(data, n, 2, r=n, tau = tau, MLEhSSALT_Obj = resMLE.1, language = "CPP",monitoring = "continuous", B = B)
row4 <- round(as.numeric(unlist(CI_i.21)[-c(9,10,12)]),2)

CI_i.22 <- CIhSSALT(data, n, 2, r=n, tau = tau, MLEhSSALT_Obj = resMLE.1, language = "CPP",monitoring = "continuous", B = B, CImethod="percentile")
row5 <- round(as.numeric(unlist(CI_i.22)[-c(9,10,12)]),2)

CI_i.23 <- CIhSSALT(data, n, 2, r=n, tau = tau, MLEhSSALT_Obj = resMLE.1, language = "CPP",monitoring = "continuous", B = B, CImethod="bca")
row6 <- round(as.numeric(unlist(CI_i.23)[-c(9,10,12)]),2)

#alpha <- 0.99
resMLE.1 <- MLEhSSALT(data,n,2,r=n,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")

CI_i.31 <- CIhSSALT(data, n, 2, r=n, tau = tau, MLEhSSALT_Obj = resMLE.1, language = "CPP",monitoring = "continuous", B = B, alpha = 0.01)
row7 <- round(as.numeric(unlist(CI_i.31)[-c(9,10,12)]),2)

CI_i.32 <- CIhSSALT(data, n, 2, r=n, tau = tau, MLEhSSALT_Obj = resMLE.1, language = "CPP",monitoring = "continuous", B = B, CImethod="percentile", alpha = 0.01)
row8 <- round(as.numeric(unlist(CI_i.32)[-c(9,10,12)]),2)

CI_i.33 <- CIhSSALT(data, n, 2, r=n, tau = tau, MLEhSSALT_Obj = resMLE.1, language = "CPP",monitoring = "continuous", B = B, CImethod="bca", alpha = 0.01)
row9 <- round(as.numeric(unlist(CI_i.33)[-c(9,10,12)]),2)

df.ii.funcs <- rbind(row1, row2, row3, row4, row5, row6, row7, row8, row9)
df.ii.funcs <- as.data.frame(df.ii.funcs, stringsAsFactors = FALSE)
colnames(df.ii.funcs) <- colnames(df.ii)[1:9]
rownames(df.ii.funcs) <- NULL
df.ii.funcs$type_c <- type_c

#Make Difference Dataframe
df.ii.diffs <- df.ii
df.ii.diffs[1:8] <- round(abs(df.ii[1:8] - df.ii.funcs[1:8]), 2)
df.ii.diffs



#############Type II - r=30
#### Given results Dataframe
th1_l <- c(8.4,18.71,18.18,4.90,16.97,16.28,0,14.71,13.76)
th1_u <- c(44.98,65.42,54.63,48.49,68.46,66.84,55.33,136.36,91.83)
th21_l <- c(0.14,0.19,0.32,0.02,0.11,0.24,0,0.01,0.06)
th21_u <- c(1.44,1.64,2.01,1.56,1.89,2.26,1.8,2.29,2.49)
th22_l <- c(2.16,3.15,3.62,1.34,2.77,3.11,0,2.25,2.38)
th22_u <- c(10.74,12.87,14.47,11.56,14.56,16.96,13.16,21.64,22.47)
p_l <-  c(0.07,0.13,0.19,0.02,0.09,0.15,0,0.05,0.08)
p_u <- c(0.67,0.63,0.68,0.73,0.66,0.72,0.84,0.73,0.78)

alpha_c <- c(rep(0.1,3),rep(0.05,3),rep(0.01,3))
type_c <- rep(c("Asymptotic", "Percentile", "BCa"), 3)

df.iii <- data.frame(
  th1_l = th1_l,
  th1_u = th1_u,
  th21_l = th21_l,
  th21_u = th21_u,
  th22_l = th22_l,
  th22_u = th22_u,
  p_l = p_l,
  p_u = p_u,
  alpha_c = alpha_c,
  type_c = type_c
)

#alpha <- 0.90
resMLE.2 <- MLEhSSALT(data,n,2,r=r,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")

CI_i.11 <- CIhSSALT(data, n, 2, r=n, tau = tau, MLEhSSALT_Obj = resMLE.2, language = "CPP",monitoring = "continuous", B = B, alpha = 0.1)
row1 <- round(as.numeric(unlist(CI_i.11)[-c(9,10,12)]),2)

CI_i.12 <- CIhSSALT(data, n, 2, r=n, tau = tau, MLEhSSALT_Obj = resMLE.2, language = "CPP",monitoring = "continuous", B = B, CImethod="percentile", alpha = 0.1)
row2 <- round(as.numeric(unlist(CI_i.12)[-c(9,10,12)]),2)

CI_i.13 <- CIhSSALT(data, n, 2, r=n, tau = tau, MLEhSSALT_Obj = resMLE.2, language = "CPP",monitoring = "continuous", B = B, CImethod="bca", alpha = 0.1)
row3 <- round(as.numeric(unlist(CI_i.13)[-c(9,10,12)]),2)

#alpha <- 0.95
resMLE.2 <- MLEhSSALT(data,n,2,r=r,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")

CI_i.21 <- CIhSSALT(data, n, 2, r=n, tau = tau, MLEhSSALT_Obj = resMLE.2, language = "CPP",monitoring = "continuous", B = B)
row4 <- round(as.numeric(unlist(CI_i.21)[-c(9,10,12)]),2)

CI_i.22 <- CIhSSALT(data, n, 2, r=n, tau = tau, MLEhSSALT_Obj = resMLE.2, language = "CPP",monitoring = "continuous", B = B, CImethod="percentile")
row5 <- round(as.numeric(unlist(CI_i.22)[-c(9,10,12)]),2)

CI_i.23 <- CIhSSALT(data, n, 2, r=n, tau = tau, MLEhSSALT_Obj = resMLE.2, language = "CPP",monitoring = "continuous", B = B, CImethod="bca")
row6 <- round(as.numeric(unlist(CI_i.23)[-c(9,10,12)]),2)

#alpha <- 0.99
resMLE.2 <- MLEhSSALT(data,n,2,r=r,tau=tau,theta21 = theta21,theta22 = theta22,p=p,language = "CPP", monitoring="continuous")

CI_i.31 <- CIhSSALT(data, n, 2, r=n, tau = tau, MLEhSSALT_Obj = resMLE.2, language = "CPP",monitoring = "continuous", B = B, alpha = 0.01)
row7 <- round(as.numeric(unlist(CI_i.31)[-c(9,10,12)]),2)

CI_i.32 <- CIhSSALT(data, n, 2, r=n, tau = tau, MLEhSSALT_Obj = resMLE.2, language = "CPP",monitoring = "continuous", B = B, CImethod="percentile", alpha = 0.01)
row8 <- round(as.numeric(unlist(CI_i.32)[-c(9,10,12)]),2)

CI_i.33 <- CIhSSALT(data, n, 2, r=n, tau = tau, MLEhSSALT_Obj = resMLE.2, language = "CPP",monitoring = "continuous", B = B, CImethod="bca", alpha = 0.01)
row9 <- round(as.numeric(unlist(CI_i.33)[-c(9,10,12)]),2)

df.iii.funcs <- rbind(row1, row2, row3, row4, row5, row6, row7, row8, row9)
df.iii.funcs <- as.data.frame(df.iii.funcs, stringsAsFactors = FALSE)
colnames(df.iii.funcs) <- colnames(df.iii)[1:9]
rownames(df.iii.funcs) <- NULL
df.iii.funcs$type_c <- type_c

#Make Difference Dataframe
df.iii.diffs <- df.iii
df.iii.diffs[1:8] <- round(abs(df.iii[1:8] - df.iii.funcs[1:8]), 2)
df.iii.diffs


# write.csv(df.i, "dfi.csv", row.names = FALSE)
# write.csv(df.i.funcs, "dfifuncs.csv", row.names = FALSE)
# write.csv(df.i.diffs, "dfidifffs.csv", row.names = FALSE)
# 
# write.csv(df.ii, "dfii.csv", row.names = FALSE)
# write.csv(df.ii.funcs, "dfiifuncs.csv", row.names = FALSE)
# write.csv(df.ii.diffs, "dfiidiffs.csv", row.names = FALSE)
# 
# write.csv(df.iii, "dfiii.csv", row.names = FALSE)
# write.csv(df.iii.funcs, "dfiiifuncs.csv", row.names = FALSE)
# write.csv(df.iii.diffs, "dfiiidiffs.csv", row.names = FALSE)