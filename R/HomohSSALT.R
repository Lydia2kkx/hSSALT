#' Homogeneity test for hSSALT
#'
#' Perform a homogeneity test under the second stress level \code{s2} of a simple hSSALT model with exponential (continuous) distribution.
#'
#' @usage HomohSSALT(data, n, censoring=1, tau, r, alpha = 0.05, M = 10000)
#'
#' @param data sample, a vector. The given data should be a censored vector with observations less than or equal to \code{n}. When censoring type is \code{2}, the length of \code{data} should be \code{r}.
#' @param n sample size, a positive integer.
#' @param censoring \code{1} for Type-I censoring or \code{2} for Type-II censoring. Default value is \code{1}.
#' @param tau If censoring type is \code{1}, \code{tau} is a vector with length 2; if censoring type is \code{2}, \code{tau} is a positive numeric value.
#' @param r If censoring type is \code{2}, \code{r} provides the pre-specified number of failures, a positive integer.
#' @param alpha significance level. Default value is \code{0.05}.
#' @param M number of simulations used to generate critical values, a positive integer. Default value is \code{10000}.
#'
#' @return An \code{hSSALTtest} object containing a hypothesis test table that reports the test statistic, the simulated critical value at the given significance level, the alternative hypothesis, and the test decision.
#'
#' @examples
#' test <- HomohSSALT(data = hSSALTdata$data, n = 35, tau = c(8, 20))
#'
#' @export



HomohSSALT <- function(data, n, censoring = 1, tau, r = NULL, alpha = 0.05, M = 10000){
  ### Part 1: Check Validity of Given Input
  ###Check the input of parameter n
  if (missing(n)) {
    stop("Error: Missing 'n'")
  }
  if (length(n)>1 || n < 0 || !is.numeric(n)) { # Added vector check
    stop("Invalid argument 'n'! The value of 'n' should be a positive number")
  }
  ###Check the input of parameter censoring
  if (!(censoring %in% c(1,2))) {
    censoring <- 1
    warning("Invalid argument 'censoring'! censoring is either 1 or 2. censoring = 1 is used instead")
  }
  ###Check the input of parameter tau
  if (missing(tau)) {
    stop("Error: Missing 'tau'")
  }
  if (!is.numeric(tau) || any(tau < 0)) { # Added !is.numeric check
    stop("Invalid argument 'tau'! The values of 'tau' should be positive numbers")
  }
  ###Check the relation between s and tau
  s <- (censoring == 1) * 2 + (censoring == 2) * 1
  if(s != length(tau)){
    stop(paste("In type", censoring, "censoring, the number of stress levels (2) and the length of tau don't match"))
  }
  ###Check the input of parameter r
  if (((censoring == 1) + !is.null(r)) == 2) {
    stop("Error: please check which censoring you want to use. If censoring is 1, r is not needed.
         If r is defined, censoring should be 2")
  }
  if (censoring == 2) {
    if(is.null(r)){
      stop("Error: Missing r")
    }else{
      if(r < 0 || !is.numeric(r) || length(r) != 1){
        stop("Invalid argument 'r'! The value of 'r' should be a positive number")
      }
    }
  }
  
  if (any(alpha < 0, alpha > 1)) {
    stop("Invalid argument 'alpha'! The value of 'alpha' should be a in range (0,1)")
  }
  
  if(any(is.integer(M), M < 0)){
    stop("Invalid argument 'M'! The value of 'M' should be a positive integer")
  }
  ############################################################################
  ############# Part 1 first stress level
  ############################################################################
  
  T1 <- data[data < tau[1]]
  n1 <- length(T1)
  T2 <- data[data >= tau[1]]
  
  if (n1 < 1){
    warning("No observation under the first stress level!")
  }
  ############################################################################
  ############# Part 2 second stress level
  ############################################################################
  
  # retrieves same result for censored and uncensored data
  n_c <- length(data)
  if(censoring == 1){
    n2 <- length(T2[T2 <= tau[2]])
    if (n > n_c){  #censored data
      t22 <- c(T2, rep(tau[2], n - n_c))
      d <- 1*(t22 < tau[2])
    }else{       #uncensored data
      d <- 1*(T2 < tau[2])              # r = tau2
      t22 <- T2 # observed or censored data
      t22[which(t22 >= tau[2])] = tau[2] #set data after censoring tau2 to tau 2
    }
  }
  
  if (censoring ==2){
    n2 <- r - n1
    cs <- T2[r - n1]
    if (n > n_c){
      t22 <- c(T2, rep(cs, n - n_c))
      d <- c(rep(1, n_c - n1), rep(0, n - n_c))
    }else{
      d <- 1*(T2 <= cs)              # censoring indicator
      t22 <- apply(cbind(T2, cs), 1, min) # observed or censored data
    }
  }
  
  n_c_data <- sum(d)
  n2_num <- n - n1
  
  ########Compute the test statistic
  obs <- sort(t22 - tau[1])
  Hs <- NULL
  for (k in 1:(n_c_data - 1)) {
    H <- -k*log(k) - (n_c_data - k)*log(n_c_data - k) + k*log(sum(obs[1:k])) + 
      (n_c_data - k)*(log(sum(obs[(k+1):n2_num])))
    Hs <- c(Hs, H)
  }
  TS <- - n_c_data*log(n_c_data) + n_c_data*log(sum(obs)) - min(Hs)
  
  ########Simulate the critical value
  TSs <- rep(0, M)
  seed <- 1
  j <- 1
  CVs <- rep(0, M)
  while (j <= M) {
    #Avner: Commented out the set.seed for now
    #set.seed(seed)
    simu_sample <- sort(rexp(n2_num))
    if(censoring == 1){
      censored_point <- -log(1 - n_c_data/n2_num)
    }else{
      censored_point <- simu_sample[r-n1]
    }
    simu_sample_censored <- apply(cbind(simu_sample, censored_point), 1, min)
    d <- c(rep(1, sum(simu_sample <= censored_point)), rep(0, (n2_num - sum(simu_sample <= censored_point))))
    n_c_sample <- sum(d)
    if(n_c_sample <= 1){
      seed <- seed + 1
      next
    }else{
      H1s <- NULL
      for (k in 1:(n_c_sample - 1)) {
        H1 <- -k*log(k) - (n_c_sample-k)*log(n_c_sample-k) + k*log(sum(simu_sample_censored[1:k])) + 
          (n_c_sample - k)*(log(sum(simu_sample_censored[(k+1):n2_num])))
        H1s <- c(H1s, H1)
      }
      TS1 <- - n_c_sample*log(n_c_sample) + n_c_sample*log(sum(simu_sample_censored)) - min(H1s)
      CVs[j] <- TS1
      seed <- seed + 1
      j <- j + 1
    }
  }
  CV_data <- quantile(CVs, probs = 1 - alpha)[[1]]
  
  ##########Decision
  decision <- if (TS > CV_data) {
    "Reject H0, data under the second stress level is heterogeneous"
  } else {
    "Cannot reject H0, data under the second stress level is homogeneous"
  }
  
  result <- list(
    statistic = c("Test Statistic" = TS),
    parameter = c("Critical Value" = CV_data, "alpha" = alpha),
    alternative = "Data are heterogeneous under second stress level",
    method = "Homogeneity Test for hSSALT",
    data.name = "Second stress level data",
    decision = decision
  )
  
  class(result) <- "hSSALTtest"
  return(result)
}

