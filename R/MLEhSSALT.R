#' Maximum Likelihood Estimation for hSSALT
#'
#' Provide point estimation of a simple hSSALT model with exponential (continuous) or geometric (interval) distribution.
#'
#' @usage MLEhSSALT(data, n, censoring = 1, tau, r = NULL, monitoring = "continuous", 
#'           delta = NULL, theta21, theta22, p, maxit = 1000, 
#'           tol = 1e-8, language = "CPP", parallel = FALSE, ncores)
#'
#' @param data sample, a vector. The given data should be a censored vector with observations less than or equal to \code{n}. When censoring type is \code{2}, the length of \code{data} should be \code{r}.
#' @param n sample size, a positive integer.
#' @param censoring \code{1} for Type-I censoring or \code{2} for Type-II censoring. Default value is \code{1}.
#' @param tau If censoring type is \code{1}, \code{tau} is a vector with length 2; if censoring type is \code{2}, \code{tau} is a positive numeric value.
#' @param r If censoring type is \code{2}, \code{r} provides the pre-specified number of failures, a positive integer.
#' @param monitoring \code{"continuous"} or \code{"interval"}. Default value is \code{"continuous"}. For interval monitoring, only equally spaced inspection is supported.
#' @param delta if interval monitoring, interval length, a positive numeric value. Default value is \code{NULL}.
#' @param theta21 initial value of \code{theta21} for the EM algorithm, can be both a numeric value or a vector of values. For an initial-value vector, the (ultimate) value with the largest log-likelihood is returned as the MLE.
#' @param theta22 initial value of \code{theta22} for the EM algorithm, can be both a numeric value or a vector of values.
#' @param p initial value of mixture proportion \code{p}, can be both a numeric value or a vector of values.
#' @param maxit The maximum number of iterations allowed, an integer. Default value is \code{1000}.
#' @param tol Tolerance limit for declaring algorithm convergence based on the change between two consecutive iterations. Default value is \code{1e-8}.
#' @param language \code{"R"} or \code{"CPP"}. Only for bootstrap methods. Default value is \code{"CPP"}.
#' @param parallel support parallel computation for multiple initial values, a logical value. Default value is \code{FALSE}.
#' @param ncores the number of cores that are used in parallelization, a positive integer.
#'
#' @seealso \code{\link{CIhSSALT}}
#'
#' @examples
#' mle <- MLEhSSALT(data = hSSALTdata$data, n = 35, censoring = 1, tau = c(8, 20),
#'        theta21 = 1, theta22 = 8, p = 0.4)
#'
#' @export


MLEhSSALT <- function( data, n, censoring = 1, tau, r = NULL, monitoring = "continuous",
                       delta = NULL, theta21, theta22, p, maxit = 1000, tol = 1e-8, 
                       language = "CPP", parallel = FALSE, ncores) {
  
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
  if ((monitoring == "interval")){ 
    
    if(censoring == 2){
      stop("Error: interval monitoring is only valid for type 1 censoring.")
    }
    ###Check the input of parameter delta
    if (is.null(delta)) {
      stop("Error: Missing 'delta'")
    }
    else if(delta < 0 || !is.numeric(delta)) { # removed ==2
      stop("Invalid argument 'delta'! The value of 'delta' should be a positive number")
    }
    if (any((tau/delta) %% 1 != 0)) {
      stop("Invalid argument 'delta'! Only equally-spaced inspection is supported. 
           The value of 'delta' should be divisible by tau")
    }
    if (!(is.integer(data))){
      stop("Error: interval monitoring is only valid for count data.")
    }
    first_interval <- as.numeric(strsplit(substring(names(data)[1],2,nchar(names(data)[1])-1), ",")[[1]])
    if (first_interval[2]-first_interval[1] != delta) {
      stop("Error: 'delta' should equal the length of the intervals in the sample.")
    }
  }
  
  ###Check the input of parameter theta21
  if (missing(theta21)) {
    stop("Error: Missing theta21")
  }
  if (any(theta21 < 0)) {
    stop("Invalid argument 'theta21'! The value of 'theta21' should be a positive vector")
  }

  ###Check the input of parameter theta22
  if (missing(theta22)) {
    stop("Error: Missing theta22")
  }
  if (any(theta22 < 0)) {
    stop("Invalid argument 'theta22'! The value of 'theta22' should be a positive vector")
  }
  
  if (missing(p)) {
    stop("Error: Missing p")
  }
  if (any(p < 0, p > 1)) {
    stop("Invalid argument 'p'! The value of 'p' should be a positive vector with entries <= 1")
  }
  
  if (!is.numeric(theta21) || !is.numeric(theta22) || !is.numeric(p)) { # Added check
    stop("Argmunts 'p', 'theta21' and 'theta22' must be positive numbers")
  }
  
  if (!(length(p)==length(theta21) && length(theta21)==length(theta22))){
    stop("Argmunts 'p', 'theta21' and 'theta22' must be of the same length")
  }
  
  if(any(c(length(p),length(theta21),length(theta22)) >1)){
    warning("Vector for initial values detected, the entries which return the largest log-likelihood are selected as the MLE")
  }
  
  if (length(tol) != 1 || length(maxit) != 1) { # Added
    stop("Arguments 'tol' and 'maxit' must not be given as vectors")
  }
  
  if (tol < 0 || !is.numeric(tol)) {
    stop("Invalid argument 'tol'! The value of 'tol' should be a positive number")
  }
  
  
  if (maxit < 0 || !is.numeric(maxit)) {
    stop("Invalid argument 'maxit'! The value of 'maxit' should be a positive number")
  }
  
  if (ncores < 0 || !is.numeric(ncores)) {
    stop("Invalid argument 'ncores'! The value of 'ncores' should be a positive integer")
  }
  if (ncores > parallel::detectCores()){
    ncores = parallel::detectCores()
    warning(paste0("System has less cores. 'ncores' is set to ", ncores))
  }
  
  if (!is.logical(parallel)){
    stop("Invalid argument 'parallel'! The value of 'parallel' should be a boolean")
  }
  
  if (!is.logical(parallel)){
    stop("Invalid argument 'parallel'! The value of 'parallel' should be a boolean")
  }
  
  if (!language %in% c("CPP", "R")){
    language <- "CPP"
    warning("Invalid argument 'language'! language is 'CPP' or 'R'. 
            language = 'CPP' is used instead")
  }
  
  # Main function
  
  if( monitoring == "continuous" ) {
    return(MLE_Exp(data = data, n = n, censoring, tau = tau, r, theta21 = theta21, 
                   theta22 = theta22, p = p, maxit, tol, language, parallel, ncores))
  } else {
    return(MLE_Geo(data = data, n = n, tau = tau, delta = delta, theta21 = theta21, 
                   theta22 = theta22, p = p, maxit, tol, language, parallel, ncores))
  }
}
