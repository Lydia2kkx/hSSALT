

CIhSSALT <- function( data , n , censoring = 1 , tau , r=NULL, monitoring = "continuous" ,
                          delta = NULL, CImethod = "asymptotic" , alpha = 0.05 , B = 1000 , theta1,
                          theta21 , theta22 , p , maxit=1000, tol=1e-8, language ="CPP", parallel=FALSE, ncores=2){
  
  ### Part 1: Check Validity of Given Input
  ###Check the input of parameter n
  if (missing(n)) {
    stop("Error: Missing 'n'")
  }
  if (n < 0 || !is.numeric(n)) {
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
  if (any(tau < 0) ) {
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
    else if((delta < 0 || !is.numeric(delta)) == 2) {
      stop("Invalid argument 'delta'! The value of 'delta' should be a positive number")
    }
    if (any((tau/delta) %% 1 != 0)) {
      stop("Invalid argument 'delta'! Only equally-spaced inspection is supported. 
           The value of 'delta' should be divisible by tau")
    }
    if (!(is.integer(data))){
      stop("Error: interval monitoring is only valid for count data.")
    }
  }
  
  ###Check the input of parameter theta21
  if (missing(theta21)) {
    stop("Error: Missing theta21")
  }
  if (any(theta21 < 0)) {
    stop("Invalid argument 'theta21'! The value of 'theta21' should be a positive vector")
  }
  # if (length(theta21) > 1) {
  #   warning("")
  # }
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
  
  if (!(length(p)==length(theta21) && length(theta21)==length(theta22))){
    stop("Argmunts 'p', 'thet21' and 'theta22' must be of the same length")
  }
  
  if(any(c(length(p),length(theta21),length(theta22)) >1)){
    warning("Vector for initial values detected, the entries which return the largest log-likelihood are selected as the MLE")
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
  
  if (any(alpha < 0, alpha > 1)) {
    stop("Invalid argument 'alpha'! The value of 'alpha' should be a in range (0,1)")
  }
  
  
  
  ###Check the input of parameter CImethod
  if(!(CImethod=="asymptotic" || CImethod=="percentile" || CImethod=="bca") ){
    CImethod <- "asymptotic"
    warning ( "’CImethod’ should be ’asymptotic’ , ’percentile’ or ’bca’!
               The default value of’ asymptotic’ is used instead.")
  }
  ###3 methods
  if(CImethod == "asymptotic"){
    CIsay_hSSALT( data, n, censoring, tau , r, monitoring, delta, alpha, theta1, p, theta21, theta22  )
  }else{
    if(CImethod == "percentile"){
      CIbs_hSSALT( data , n , censoring , tau, r , monitoring , delta , alpha , B,
                   theta21 , theta22 , p , maxit , tol , language , parallel , ncores )
    }else{
      CIbca_hSSALT( data , n , censoring , tau , r, monitoring , alpha , B,
                    theta21 , theta22 , p , maxit , tol , language , parallel , ncores )
    }
  }
}