source("rhSSALT.R")
set.seed(1)
n <- 40
sample <- rhSSALT(n, 1, tau = c(8, 16), theta1 = 33.12, theta21 = 2, theta22 = 15,
                  p = 0.4, monitoring = "interval", delta = 1)
sample$Censored_dat

jackknife_samples_interval <- list()
index <- 1
for (bin in seq_along(sample$Censored_dat)) {
  if (sample$Censored_dat[bin] > 0) {
    for (k in 1:sample$Censored_dat[bin]) {
      temp <- sample$Censored_dat
      temp[bin] <- temp[bin] - 1
      jackknife_samples_interval[[index]] <- temp
      index <- index + 1
    }
  }
}

jackknife_samples_interval
