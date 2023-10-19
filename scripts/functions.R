# Regression diagnostics for mids object
mids_diagnostics <- function(model, name = "") {
  res <- sapply(model$analyses, residuals)
  res <- apply(res, 1, mean)
  fit <- sapply(model$analyses, fitted)
  fit <- apply(fit, 1, mean)
  res_standard <- sapply(model$analyses, rstandard)
  res_standard <- apply(res_standard, 1, mean)
  
  plot(fit, res, main = paste(name))
  qqnorm(res_standard, main = paste(name))
  qqline(res_standard)
}
