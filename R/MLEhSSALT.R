MLEhSSALT <- function( data , n , censoring=1 , tau , r=NULL , monitoring="continuous" ,
                               delta=NULL , theta21 , theta22 , p, maxit=1000, tol=1e-8, language ="CPP", parallel=FALSE, ncores=2) {
  
  
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
  
  
  
  
  # if (length(theta22) > 1) {
  #   theta22 <- theta22[1]
  #   warning("Invalid argument 'theta22'! 'theta22' should a number. Only the first element in the vector is used instead")
  # }
  
  ### 
  # if( !is.vector( data ) ) {
  #   stop( " Invalid argument ’data’! The type of ’data’ should be a vector " )
  # }
  # if( any(data < 0) || !is.numeric (data) ) {
  #   stop( " Invalid argument ’data’ ! The values of ’data’ should be positive numbers " )
  # }
  
  
  
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
    return(MLE_Exp( data , n , censoring , tau , r , theta21 , theta22 , p, maxit, tol, language, parallel, ncores ))
  } else {
    return(MLE_Geo( data=data , n=n , tau=tau , delta=delta , theta21=theta21 , theta22=theta22 , p=p, maxit, tol, language, parallel, ncores ))
  }
}



# 
# ###Check the input of parameter n
# if( n < 0 || !is.numeric(n) ) {
#   stop( "Invalid argument ’n’! The value of ’n’ should be a postive number" )
# }
# if ( n < length(data) ) {
#   stop( " Invalid argument 'n' ! The value of 'n' should be larger  than the
#            length of data " )
# }
# ###Check the input of parameter censoring
# if((censoring == 1) + (censoring == 2) < 1 ) {
#   censoring <- 1
#   warning ( " Invalid argument ’ censoring ’ ! censoring is either 1 or 2 .
#            censoring = 1 is used instead" )
# }
# 
# ###Check the input of parameter monitoring
# if((monitoring == "continuous") + (monitoring == "interval") < 1 ) {
#   censoring <- "continuous"
#   warning ( " Invalid argument ’ censoring ’ ! censoring is either 'continuous' or 'interval' .
#            continuous monitoring is used instead" )
# }
# 
# # Check for additional parameters
# #
# # d
# if(is.element(FALSE,unique(d)==c(1,0))){
#   stop(" Invalid argument 'd' ! The entries of 'd' should either be zero or 1" )
# }
# 
# # if(length(d)!=n){
# #   stop(" Invalid argument 'd' ! d must have n entries" )
# # }
# 
# ###Check the input of parameter tau
# if( tau < 0 || !is.numeric(tau) ) {
#   stop( "Invalid argument ’tau’! The value of ’tau’ should be a postive number" )
# }
# 
# ###Check the input of parameter r
# if( r < 0 || !is.numeric(r) ) {
#   stop( "Invalid argument ’r’! The value of ’r’ should be a postive number" )
# }
# 
# # if(censoring==2 && !is.integer(r)) {
# #   stop( "Invalid argument ’r’! In type-II censoring The value of ’r’ should be an integer" )
# # }
# 
# ###Check the input of parameter delta
# if( delta < 0 || !is.numeric(delta) || delta > 1 ) {
#   stop( "Invalid argument ’delta’! The value of ’delta’ should be a postive number in (0,1)" )
# }
# 
# 
# ###Check the input of parameter p
# if(is.element(TRUE, p <0 ) || is.element(TRUE, p > 1) ) {
#   stop( "Invalid argument ’p’! The value of ’p’ should be a postive number in (0,1)" )
#   
# }
# 
# 
# ###Check the input of parameter theta21
# if( is.element(FALSE, theta21 > 0) ) {
#   stop( "Invalid argument ’theta21’! The value of ’theta21’ should be  postive" )
# }
# 
# ###Check the input of parameter theta21
# if( is.element(FALSE, theta22 > 0) ) {
#   stop( "Invalid argument ’theta22’! The value of ’theta22’ should be  postive" )
# }
# 



