#' @export
print.hSSALTMLE <- function(x, ...) {

  if (x$message == "convergent") {
    cat(sprintf("The EM algorithm converged after %s iterations.\n", x$iteration))
  } else {
    cat("The EM algorithm failed to converge.\n")
  }
  
  cat("Log-likelihood:",
      formatC(x$loglik, digits = 6, format = "f"), "\n\n")
  
  cat("Maximum Likelihood Estimates:\n")
  ests <- formatC(unlist(x$mle), digits = 6, format = "f")
  print(ests, quote = FALSE)
  cat("\n")
  
  # To be able to capture list in variable with print
  invisible(x)
}