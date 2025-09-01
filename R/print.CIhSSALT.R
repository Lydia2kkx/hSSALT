#' @export
print.CIhSSALT <- function(x, ...) {
  
  title <- sprintf("%s confidence intervals for theta1, theta21, theta22 and p", x$type)
  
  # Help functions
  number_formatting <- function(z, dig=6) formatC(z, digits = dig, format = "f", drop0trailing = FALSE)
  separation_line <- paste(rep("-", 50), collapse = "")
  
  # Print
  cat(title, "\n", sep = "")
  cat(sprintf("Significance level: alpha = %s\n", number_formatting(x$alpha, dig=2)))
  
  if (x$type != "Asymptotic") {
    cat(sprintf("Valid bootstrap replications: B = %s\n", x$B))
    cat(sprintf("Total bootstrap replications: j = %s\n", x$j)) 
  }
  
  cat("\n", separation_line, "\n", sep = "")
  
  cat(sprintf("theta1  (%s, %s)\n",  number_formatting(x$theta1[1]),  number_formatting(x$theta1[2])))
  cat(sprintf("theta21 (%s, %s)\n", number_formatting(x$theta21[1]), number_formatting(x$theta21[2])))
  cat(sprintf("theta22 (%s, %s)\n", number_formatting(x$theta22[1]), number_formatting(x$theta22[2])))
  cat(sprintf("p       (%s, %s)\n", number_formatting(x$p[1]), number_formatting(x$p[2])))
  
  # To be able to capture list in variable with print
  invisible(x)
}