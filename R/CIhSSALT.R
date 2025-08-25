#' Confidence Intervals for hSSALT
#'
#' Provide interval estimation of a simple hSSALT model with exponential (continuous) or geometric (interval) distribution.
#'
#' @param data sample, a vector. The given data should be a censored vector with observations less than or equal to \code{n}. When censoring type is \code{2}, the length of \code{data} should be \code{r}.
#' @param n sample size, a positive integer.
#' @param MLEhSSALT an \code{MLEhSSALT} object, returned by \code{MLEhSSALT()}.
#' @param censoring \code{1} for Type-I censoring or \code{2} for Type-II censoring. Default value is \code{1}.
#' @param tau If censoring type is \code{1}, \code{tau} is a vector with length 2; if censoring type is \code{2}, \code{tau} is a positive numeric value.
#' @param r If censoring type is \code{2}, \code{r} provides the pre-specified number of failures, a positive integer.
#' @param monitoring \code{"continuous"} or \code{"interval"}. Default value is \code{"continuous"}. For interval monitoring, only equally spaced inspection is supported.
#' @param delta if interval monitoring, interval length, a positive numeric value. Default value is \code{NULL}.
#' @param CImethod \code{"asymptotic"}, \code{"percentile"} or \code{"bca"} for asymptotic CIs, bootstrap percentile CIs and bootstrap bias-corrected and accelerated (BCa) bootstrap intervals. Default value is \code{"asymptotic"}.
#' @param alpha significance level. Default value is \code{0.05}.
#' @param B number of bootstrap repetitions, a positive integer, default value is \code{1000}.
#' @param maxit The maximum number of iterations allowed, a positive integer. Only for bootstrap methods. Default value is \code{1000}.
#' @param tol Tolerance limit for declaring algorithm convergence based on the change between two consecutive iterations. Only for bootstrap methods. Default value is \code{1e-8}.
#' @param language \code{"R"} or \code{"CPP"}. Only for bootstrap methods. Default value is \code{"CPP"}.
#' @param parallel support parallel computation, a logical value. Only for bootstrap methods. Default value is \code{FALSE}.
#' @param ncores the number of cores that are used in parallelization, a positive integer.
#'
#' @return A \code{CIhSSALT} object that includes the type of returned CIs and the CIs for four parameters at a given significance level.
#'
#' @examples
#' sample <- rhSSALT(n = 30, tau = c(5, 10), theta1 = 10, theta21 = 5, theta22 = 8, p = 0.4)
#' MLE <- MLEhSSALT(data = sample$`censored sample`, n = 30, censoring = 1, tau = c(5, 10), theta21 = 5, theta22 = 8, p = 0.4)
#' ci <- CIhSSALT(data = sample$`censored sample`, n = 30, MLEhSSALT = MLE, tau = c(5, 10))
#'
#' @export



CIhSSALT <- function(data, n, MLEhSSALT_Obj, censoring = 1, tau, r=NULL, monitoring = "continuous",
                     delta = NULL, CImethod = "asymptotic", alpha = 0.05, B = 1000, maxit = 1000, 
                     tol=1e-8, language = "CPP", parallel = FALSE, ncores = 2, grid = FALSE){
  
  theta1 <- MLEhSSALT_Obj$mle$theta1
  theta21 <- MLEhSSALT_Obj$mle$theta21
  theta22 <- MLEhSSALT_Obj$mle$theta22
  p <- MLEhSSALT_Obj$mle$p1
  
  ### Part 1: Check Validity of Given Input
  ###Check the input of parameter n
  if (missing(n)) {
    stop("Error: Missing 'n'")
  }
  if (length(n)>1 || n < 0 || !is.numeric(n)) {
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
  if (!is.numeric(tau) || any(tau < 0) ) {
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
    else if(delta < 0 || !is.numeric(delta)) {
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
  
  if (any(alpha < 0, alpha > 1)) {
    stop("Invalid argument 'alpha'! The value of 'alpha' should be a in range (0,1)")
  }
  
  
  
  ###Check the input of parameter CImethod
  if(!(CImethod=="asymptotic" || CImethod=="percentile" || CImethod=="bca") ){
    CImethod <- "asymptotic"
    warning ( "'CImethod' should be 'asymptotic' , 'percentile' or 'bca'!
               The default value of 'asymptotic' is used instead.")
  }
  
  ###Check MLE Object
  if (class(MLEhSSALT_Obj) != "hSSALTMLE") {
    stop("Expected an object of class 'hSSALTMLE'.")
  }
  if (MLEhSSALT_Obj$message != "convergent") {
    warning ("The EM Algorithm failed to converge! Values of 'p', 'theta21', 'theta22' may not equal the MLEs.")
  }
  if (is.na(MLEhSSALT_Obj$mle$theta1)) {
    warning("No value for 'theta1' detected.")
  } else {
    if (MLEhSSALT_Obj$mle$theta1 < MLEhSSALT_Obj$mle$theta22) {
      warning("Mean lifetime in first stress level must be greater than mean lifetimes in the second stress level.")
    }
  }
  
  ###3 methods
  if(CImethod == "asymptotic"){
    CIsay_hSSALT(data, n, censoring, tau, r, monitoring, delta, alpha, theta1, p, theta21, theta22)
  }else{
    if(CImethod == "percentile"){
      CIbs_hSSALT(data, n, censoring, tau, r, monitoring, delta, alpha, B, theta1,
                  theta21, theta22, p, maxit, tol, language, parallel, ncores, grid)
    }else{
      CIbca_hSSALT(data, n, censoring, tau, r, monitoring, delta, alpha, B,
                   theta1, theta21, theta22, p, maxit, tol, language, parallel, ncores, grid)
    }
  }
}
