#' @export
print.hSSALTtest <- function(x, ...) {
  
  cat("\n", x$method, "\n\n", sep = "")
  cat("data:  ", x$data.name, "\n", sep = "")
  cat(sprintf("%s = %.4f", names(x$statistic), x$statistic), ", ",
      paste(names(x$parameter), "=", format(x$parameter, digits = 5), collapse = ", "),
      "\n", sep = "")
  cat("alternative hypothesis: ", x$alternative, "\n", sep = "")
  
  ###Decision
  cat("\nDecision: ", x$decision, "\n", sep = "")
}