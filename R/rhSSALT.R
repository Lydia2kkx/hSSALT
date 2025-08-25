#' Simulate a simple hSSALT random dataset
#'
#' Simulate a simple hSSALT random dataset with exponential (continuous) or geometric (interval) distribution.
#'
#' @param n sample size, an integer.
#' @param censoring \code{1} for Type-I censoring or \code{2} for Type-II censoring. Default value is \code{1}.
#' @param tau If censoring type is \code{1}, \code{tau} is a vector with length 2; if censoring type is \code{2}, \code{tau} is a positive numeric value.
#' @param r If censoring type is \code{2}, \code{r} provides the pre-specified number of failures, a positive integer.
#' @param monitoring \code{continuous} or \code{interval}, default value is \code{continuous}. For interval monitoring, only equally spaced inspection is supported.
#' @param delta if interval monitoring, interval length, a positive numeric value. Default value is \code{NULL}.
#' @param theta1 mean lifetime parameter in the exponential distribution under \code{s1}, a numeric value.
#' @param theta21 mean lifetime parameter in the exponential distribution of the first group under \code{s2}, a positive numeric value.
#' @param theta22 mean lifetime parameter in the exponential distribution of the second group under \code{s2}, a positive numeric value.
#' @param p mixture proportion, a numeric value between 0 and 1.
#'
#' @return A list consisting of four sub-lists: censored sample, the observed number of censored failures under \code{s1} and \code{s2}, complete sample, the observed number of failures under \code{s1} and \code{s2}.
#'
#' @examples
#' sample <- rhSSALT(n = 30, tau = c(5, 10), theta1 = 10, theta21 = 5, theta22 = 8, p = 0.4)
#'
#' @export

rhSSALT <- function(n, censoring = 1, tau, r = NULL, monitoring = "continuous", delta = NULL, 
                    theta1, theta21, theta22, p){
  ### Part 1: Check Validity of Given Input
  ###Check the input of parameter n
  if (missing(n)) {
    stop("Error: Missing 'n'")
  }
  if (n < 0 || !is.numeric(n)) {
    stop("Invalid argument 'n'! The value of 'n' should be a positive number")
  }
  ###Check the input of parameter censoring
  if ((censoring == 1) + (censoring == 2) < 1) {
    censoring <- 1
    warning("Invalid argument 'censoring'! censoring is either 1 or 2. censoring = 1 is used instead")
  }
  ###Check the input of parameter tau
  if (missing(tau)) {
    stop("Error: Missing 'tau'")
  }
  if (any(tau < 0, !is.numeric(tau))) {
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
  ###Check the input of parameter monitoring
  if (!monitoring %in% c("continuous", "interval")){
    monitoring <- "continuous"
    warning("Invalid argument 'monitoring'! monitoring is 'continuous' or 'interval'. 
            monitoring = 'continuous' is used instead")
  }
  ######The following check is newly added.
  if ((monitoring == "continuous") + (!is.null(delta))  == 2){
    warning("Conflict in argument 'monitoring' and delta'! monitoring = 'continuous' is used instead")
  }
  if ((monitoring == "interval") + (censoring == 2) == 2){
    stop("Error: interval monitoring is only valid for type 1 censoring.")
  }
  ###Check the input of parameter delta
  if ((monitoring == "interval") + (is.null(delta)) == 2) {
    stop("Error: Missing 'delta'")
  }
  if ((monitoring == "interval") + (delta < 0 || !is.numeric(delta)) == 2) {
    stop("Invalid argument 'delta'! The value of 'delta' should be a positive number")
  }
  if (any((tau/delta) %% 1 != 0)) {
    stop("Invalid argument 'delta'! Only equally-spaced inspection is supported. 
         The value of 'delta' should be divisible by tau")
  }
  ###Check the input of parameter theta1
  if (missing(theta1)) {
    stop("Error: Missing theta1")
  }
  if (any(theta1 < 0, !is.numeric(theta1))) {
    stop("Invalid argument 'theta1'! The value of 'theta1' should be a positive number")
  }
  ###Check the input of parameter theta21
  if (missing(theta21)) {
    stop("Error: Missing theta21")
  }
  if (any(theta21 < 0, !is.numeric(theta21))) {
    stop("Invalid argument 'theta21'! The value of 'theta21' should be a positive number")
  }
  if (length(theta21) > 1) {
    theta21 <- theta21[1]
    warning("Invalid argument 'theta21'! 'theta21' should a number. Only the first element in the vector is used instead")
  }
  ###Check the input of parameter theta22
  if (missing(theta22)) {
    stop("Error: Missing theta22")
  }
  if (any(theta22 < 0, !is.numeric(theta22))) {
    stop("Invalid argument 'theta22'! The value of 'theta22' should be a positive number")
  }
  if (length(theta22) > 1) {
    theta22 <- theta22[1]
    warning("Invalid argument 'theta22'! 'theta22' should a number. Only the first element in the vector is used instead")
  }
  ###Check the input of parameter p
  if (missing(p)) {
    stop("Error: Missing p")
  }
  if (any(p < 0, p > 1, (!is.numeric(p)))) {
    stop("Invalid argument 'p'! The value of 'p' should between 0 and 1")
  }
  if (length(p) > 1) {
    p <- p[1]
    warning("Invalid argument 'p'! 'p' should a number. Only the first element in the vector is used instead")
  }
  
  ### Part 2: Simulate Data According to the Input
  T1_simu <- rexp(n, rate = 1/theta1)
  n1 <- sum(T1_simu <= tau[1])
  T1 <- sort(T1_simu[T1_simu <= tau[1]])
  n21 <- rbinom(1, (n-n1), p)
  T21 <- -theta21*log(1-runif(n21)) + tau[1]
  n22 <- n-n1-n21
  T22 <- -theta22*log(1-runif(n22)) + tau[1]
  T2 <- sort(c(T21, T22))
  Simu_T <- sort(c(T1, T2))
  if(censoring == 1){
    Simu_T_Censored <- Simu_T[Simu_T <= tau[2]]
    if(monitoring == "continuous"){
      return(list(Censored_dat = Simu_T_Censored, Censored_num_level = c(n1, length(Simu_T_Censored) - n1), 
                  Full_dat = Simu_T, Full_num_level = c(n1, n-n1)))
    }else{
      Simu_njs <- table(cut(Simu_T, breaks = seq(0, ceiling(max(Simu_T)/delta)*delta, delta)))
      Simu_njs_Censored <- table(cut(Simu_T_Censored, breaks = seq(0, tau[2], delta)))
      return(list(Censored_dat = Simu_njs_Censored, Censored_num_level = c(n1, sum(Simu_njs_Censored) - n1), 
                  Full_dat = Simu_njs, Full_num_level = c(n1, n-n1)))
    }
  }else{
    Simu_T_Censored <- Simu_T[1:r]
    return(list(Censored_dat = Simu_T_Censored, Censored_num_level = c(n1, length(Simu_T_Censored) - n1), 
                Full_dat = Simu_T, Full_num_level = c(n1, n-n1)))
  }
}

