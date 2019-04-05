#' @title Distribution of future treatment effects.
#' @author Mohamed Amine BAYAR
#' \code{DistributionFutureTreatmentEffects} derives the distribution of future treatment effects in terms of hazard rate.
#' 
#' @param index An integer specifying the index of the distribution.
#'      1 : D1 the historical distribution
#       2 : D2 the optimistic distribution
#       3 : D3 the pessimistic distribution
#       4 : D4 the very pessimistic distribution
#' @param lambda_0 A positive real number specifying the hazard rate associated with control arm of the first trial of the series.
#' 
#' @example
#' DistributionFutureTreatmentEffects(index = 1, lambda_0 = log(2))
#' 
#' @return a vector with three elements (a,b,sigma)
#' characterizing the distribution of future treatment effects in terms of hazard rate.

DistributionFutureTreatmentEffects <- function(index, lambda_0) {
  
  d <- 5 ## We considered that the average duration of a clinical trial is 5 years
  if (index==1) {
    e   <- 0.95
    p   <- 0.02
  } 
  if (index==2) {
    e   <- 0.925
    p   <- 0.02
  } 
  if (index==3) {
    e   <- 0.975
    p   <- 0.01
  } 
  if (index==4) {
    e   <- 1
    p   <- 0.01
  }  
  z        <- qnorm(p,0,1)
  aa       <- log(e)/d
  sigma    <- (log(0.5)-log(e))/(sqrt(2)*z)
  bb       <- log(e) + log(lambda_0) - 0.5*sigma^2
  return(c(round(aa,4), round(bb,3), round(sigma,3)) )
}
