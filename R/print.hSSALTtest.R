#' @export
print.hSSALTtest <- function(x, ...) {
  # Standard htest output style
  cat("\n", x$method, "\n\n", sep = "")
  cat("data:  ", x$data.name, "\n", sep = "")
  cat(sprintf("%s = %.4f", names(x$statistic), x$statistic), ", ",
      paste(names(x$parameter), "=", format(x$parameter, digits = 5), collapse = ", "),
      "\n", sep = "")
  cat("alternative hypothesis: ", x$alternative, "\n", sep = "")
  # Add decision
  cat("\nDecision: ", x$decision, "\n", sep = "")
}