#Compute the log ratio to tackle composionality
normalise_100 <- function(x) {(x/rowSums(x))*100}

log_ratio <- function(data)
{
  # Compute the logarithmus
  log_data <- log(data
                +1)
  # Calculate exponential function of column-wise mean values for finite log transformed data
  gm <- exp(mean(log_data[is.finite(log_data)]))
  # Compute the logarithmus
  log_gm <- log(gm)
  # Take the difference of both log-transformed datasets
  data <- log_data - log_gm
  # Return the new OTU table
  return(data)
}

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}