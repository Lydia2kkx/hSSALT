

#Continuous
theta1_l <- c(df.cont.asymptotic.cpp$theta1_l - df.cont.asymptotic.R$theta1_l, df.cont.percentile.cpp$theta1_l - df.cont.percentile.R$theta1_l, df.cont.bca.cpp$theta1_l - df.cont.bca.R$theta1_l)
theta1_u <- c(df.cont.asymptotic.cpp$theta1_u - df.cont.asymptotic.R$theta1_u, df.cont.percentile.cpp$theta1_u - df.cont.percentile.R$theta1_u, df.cont.bca.cpp$theta1_u - df.cont.bca.R$theta1_u)

theta21_l <- c(df.cont.asymptotic.cpp$theta21_l - df.cont.asymptotic.R$theta21_l, df.cont.percentile.cpp$theta21_l - df.cont.percentile.R$theta21_l, df.cont.bca.cpp$theta21_l - df.cont.bca.R$theta21_l)
theta21_u <- c(df.cont.asymptotic.cpp$theta21_u - df.cont.asymptotic.R$theta21_u, df.cont.percentile.cpp$theta21_u - df.cont.percentile.R$theta21_u, df.cont.bca.cpp$theta21_u - df.cont.bca.R$theta21_u)

theta22_l <- c(df.cont.asymptotic.cpp$theta22_l - df.cont.asymptotic.R$theta22_l, df.cont.percentile.cpp$theta22_l - df.cont.percentile.R$theta22_l, df.cont.percentile.cpp$theta22_l - df.cont.percentile.R$theta22_l)
theta22_u <- c(df.cont.asymptotic.cpp$theta22_u - df.cont.asymptotic.R$theta22_u, df.cont.percentile.cpp$theta22_u - df.cont.percentile.R$theta22_u, df.cont.percentile.cpp$theta22_u - df.cont.percentile.R$theta22_u)

p_l <- c(df.cont.asymptotic.cpp$p_l - df.cont.asymptotic.R$p_l, df.cont.percentile.cpp$p_l - df.cont.percentile.R$p_l, df.cont.bca.cpp$p_l - df.cont.bca.R$p_l)
p_u <- c(df.cont.asymptotic.cpp$p_u - df.cont.asymptotic.R$p_u, df.cont.percentile.cpp$p_u - df.cont.percentile.R$p_u, df.cont.bca.cpp$p_u - df.cont.bca.R$p_u)

type <- c(noquote(rep("Asymptotic", 20)), noquote(rep("Percentile", 20)), noquote(rep("BCa", 20)))
seed <- rep(1:20, 3)

cont.df <- data.frame(seed, theta1_l,theta1_u,theta21_l, theta21_u, theta22_l, theta22_u, p_l, p_u, type)
write.csv(cont.df, "testing_diffs_cont_adjusted.csv", row.names=FALSE)





#####Interval
theta1_l <- c(df.int.asymptotic.cpp$theta1_l - df.int.asymptotic.R$theta1_l, df.int.percentile.cpp$theta1_l - df.int.percentile.R$theta1_l, df.int.bca.cpp$theta1_l - df.int.bca.R$theta1_l)
theta1_u <- c(df.int.asymptotic.cpp$theta1_u - df.int.asymptotic.R$theta1_u, df.int.percentile.cpp$theta1_u - df.int.percentile.R$theta1_u, df.int.bca.cpp$theta1_u - df.int.bca.R$theta1_u)

theta21_l <- c(df.int.asymptotic.cpp$theta21_l - df.int.asymptotic.R$theta21_l, df.int.percentile.cpp$theta21_l - df.int.percentile.R$theta21_l, df.int.bca.cpp$theta21_l - df.int.bca.R$theta21_l)
theta21_u <- c(df.int.asymptotic.cpp$theta21_u - df.int.asymptotic.R$theta21_u, df.int.percentile.cpp$theta21_u - df.int.percentile.R$theta21_u, df.int.bca.cpp$theta21_u - df.int.bca.R$theta21_u)

theta22_l <- c(df.int.asymptotic.cpp$theta22_l - df.int.asymptotic.R$theta22_l, df.int.percentile.cpp$theta22_l - df.int.percentile.R$theta22_l, df.int.percentile.cpp$theta22_l - df.int.percentile.R$theta22_l)
theta22_u <- c(df.int.asymptotic.cpp$theta22_u - df.int.asymptotic.R$theta22_u, df.int.percentile.cpp$theta22_u - df.int.percentile.R$theta22_u, df.int.percentile.cpp$theta22_u - df.int.percentile.R$theta22_u)

p_l <- c(df.int.asymptotic.cpp$p_l - df.int.asymptotic.R$p_l, df.int.percentile.cpp$p_l - df.int.percentile.R$p_l, df.int.bca.cpp$p_l - df.int.bca.R$p_l)
p_u <- c(df.int.asymptotic.cpp$p_u - df.int.asymptotic.R$p_u, df.int.percentile.cpp$p_u - df.int.percentile.R$p_u, df.int.bca.cpp$p_u - df.int.bca.R$p_u)

type <- c(noquote(rep("Asymptotic", 20)), noquote(rep("Percentile", 20)), noquote(rep("BCa", 20)))
seed <- rep(1:20, 3)

intervals.df <- data.frame(seed, theta1_l,theta1_u,theta21_l, theta21_u, theta22_l, theta22_u, p_l, p_u, type)
write.csv(intervals.df, "testing_diffs_interval_adjusted.csv", row.names=FALSE)