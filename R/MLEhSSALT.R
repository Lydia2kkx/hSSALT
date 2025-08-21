MLEhSSALT <- function( data, n, censoring = 1, tau, r = NULL, monitoring = "continuous",
                       delta = NULL, theta21, theta22, p, maxit = 1000, tol = 1e-8, 
                       language = "CPP", parallel = FALSE, ncores = 2) {
  
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
