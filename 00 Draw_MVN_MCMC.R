#rm(list = ls())

library(mvtnorm)

draw.mvn.mcmc <- function(n, d = 1, mu = 0, sigma = 1, cand.s = 0.5, bias = 0)
{
  
  mu <- matrix(mu, ncol = 1) + bias
  sigma <- as.matrix(sigma)
  
  acp <- 0
  
  ret <- matrix((sample(5:10, 1) + rnorm(d)), ncol = d, nrow = n, byrow = T)
  
  for (i in 2:n){
    cand.samp <- rmvnorm(1, ret[i-1, ], cand.s^2 * diag(d))
    
    HR <- dmvnorm(cand.samp, mu, sigma) / dmvnorm(ret[i-1, ], mu, sigma)
    
    U <- runif(1)
    if (U < HR){
      ret[i, ] <- cand.samp
      acp <- acp + 1
    }else{
      ret[i, ] <- ret[i-1, ]
    }
    
  }
  
  #-------------
  
  return(list(samp = ret, acp.rat = (acp/n)))
}
