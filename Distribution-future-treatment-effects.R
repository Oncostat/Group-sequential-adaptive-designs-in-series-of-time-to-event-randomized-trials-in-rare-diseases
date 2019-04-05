#################################################################
#             Function distrib: HR distribution                 #
#                       input : Index                           #
#                 output : a list (Index, Mean, STD)            #
#################################################################

distrib     <- function(index, lambda_0)
                {
                d                 <- 5
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
