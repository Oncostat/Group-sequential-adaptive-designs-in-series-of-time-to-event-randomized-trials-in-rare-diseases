#' @title Simes test for the first stage and for the second stage.
#' @author Mohamed Amine BAYAR
#' \code{SimesFirstStage SimesSecondStage} computes the a Simes test-adjusted p-value.
#' 
#' @param x A vector of p-values.
#'
#' @return A Simes test-adjusted p-value

SimesFirstStage <- function(x){
  y <- numeric()
  x <- apply(x, 2, sort)
  for (i in 1:ncol(x))
    y[i] <- min(x[,i]*(length(x[,i])/(1:length(x[,i]))))
  return(y)
}

SimesSecondStage <- function(x){
  y <- numeric()
  x <- apply(x, 2, sort)
  for (i in 1:ncol(x))
    y[i] <- min(x[,i])
  return(y)
}

#' @title Inverse normal comination function.
#' @author Mohamed Amine BAYAR
#' \code{InverseNormalCombination} computes the combined p-value using inverse normal combination function.
#' 
#' @param p1,p2 Two p-values.
#' @param w     Weight of the first stage.
#'
#' @return A combined p-value using inverse normal combination function.
#' 
#' 
InverseNormalCombination <- function(p1,p2,w){
  qnormal1 <- qnorm(1-p1)
  qnormal2 <- qnorm(1-p2)
  pnormal  <- pnorm(sqrt(w)*qnormal1 + sqrt(1-w)*qnormal2) 
  return(1-pnormal)
}
