#rm(list = ls())

#===============================================================================
##functions for stein thinning

#calculating the pairwise-KSD
Kp.xy <- function(x, y, n, b = -0.5, c = 1)
{
  
  d <- length(x)
  w.sq <- numeric()
  
  b.x <- - x
  b.y <- - y
  
  val <-  c^2 + norm(x-y, type = "2")^2 
  
  for ( j in 1:d){
    
    term1 <- b.x[j] * b.y[j] * (  val^b  )
    term2 <- b.x[j] * ( b * val^(b-1) * (-2) * (x[j] - y[j])  )
    term3 <- b.y[j] * (  b * val^(b-1) * 2 * (x[j] - y[j])  )
    term4 <-(  - 2 * b * (  (b-1) * val^(b-2) * 2 * (x[j] - y[j])^2  + val^(b-1)  )  )
    
    w.sq[j] <-  term1 + term2 + term3 + term4
    
  }
  
  ret <- sum(w.sq)
  
  return(ret)
}

#--  --  --  --  --  --  --  --  --  --  -- 
#choosing the indices

stein.thin.indices <- function(samp, m){
  
  n <- nrow(samp)
  
  #calculating KSD of the form K_p(X_i, x_i)
  KSD.n <- numeric()
  for (i.n in 1:n)
  { 
    KSD.n[i.n] <- Kp.xy(samp[i.n,],samp[i.n,], n)
  }
  
  #We are intended to choose the indices
  ind.chosen <- numeric()
  
  ind.chosen[1] <- which.min(KSD.n)    #first value is based on KSD.n only
  
  store <- numeric()
  #loop for j = 2:m
  for (j in 2:m){
    
    KSD.trail.n <- numeric()           #storing the required quantity for each n
    
    #loop for n
    for (i in 1:n){          
      f.term <- KSD.n[i] / 2           #first term
      
      s.term <- 0                      #second term
      for (j.d in 1:(j-1))             #loop for j'
      {            
        s.term <- s.term + Kp.xy( samp[ind.chosen[j.d],], samp[i,], n)
        
      }
      
      store <- c(store, s.term)
      
      KSD.trail.n[i] <- f.term + s.term
    }
    
    ind.chosen[j] <- which.min(KSD.trail.n)
    
  }
  
  return(ind.chosen)
}


#===============================================================================
##the main function

source("00 Draw_MVN_MCMC.R")

thin.comp <- function(rep = 1e2, n.all = c(50, 100, 200, 500, 1000)){
  
  #true values
  d = 1
  mu.true <- rep(0,d)
  sigma.true <- diag(d)
  
  #storage for the values
  store.arr <- array(0, dim = c(length(n.all), 4, rep))
  
  #loops start here
  for (i.n in 1:length(n.all)){
    
    for (i.rep in 1:rep){
      
      set.seed(100 + i.n + i.rep)
      
      n.cur <- n.all[i.n]
      m <- ceiling(log(n.cur)^2)
      
      samp <- as.matrix(draw.mvn.mcmc(n.cur, d = d, mu = mu.true, sigma = sigma.true, bias = 2)[[1]])
      
      #---------------------------------------------
      #random thinning
      ind.1 <- sample(1:n.cur, m, replace = T)
      samp.1 <- as.matrix(samp[ind.1, ])
      mu.est.1 <- colMeans(samp.1)
      store.arr[i.n, 1, i.rep] <- norm(mu.est.1 - mu.true, type = "2")
      
      #---------------------------------------------
      #fixed thinning
      ind.2 <- seq(1, m * floor(n.cur/m), by = floor(n.cur/m))
      samp.2 <- as.matrix(samp[ind.2, ])
      mu.est.2 <- colMeans(samp.2)
      store.arr[i.n, 2, i.rep] <- norm(mu.est.2 - mu.true, type = "2")
      
      #---------------------------------------------
      #Stein thinning
      ind.3 <- stein.thin.indices(samp, m)
      samp.3 <- as.matrix(samp[ind.3, ]) 
      mu.est.3 <- colMeans(samp.3)
      store.arr[i.n, 3, i.rep] <- norm(mu.est.3 - mu.true, type = "2")
      
      #---------------------------------------------
      #No thinning
      mu.est.4 <- colMeans(samp)
      store.arr[i.n, 4, i.rep] <- norm(mu.est.4 - mu.true, type = "2")
      
      #---------------------------------------------
      #aesthetics
      print("-------------------------------")
      print(paste0("We are at n = ",n.all[i.n]," and rep = ", i.rep))
      
    }
    
  }
  
  return(store.arr)
  
}

#===============================================================================

out <- thin.comp(rep = 25, n.all = c(50, 100, 200, 300, 400, 500))
save(out, file = "05 Thin_comp_Biased MCMC.Rdata")

#----------------------
#rm(list = ls())

load("05 Thin_comp_Biased MCMC.Rdata")
tab <- as.table(round(apply(out, c(1, 2), mean),5))
colnames(tab) <- c("Fixed_Thinning", "Random_Thinning", "Stein_Thinning", "Full Sample")
rownames(tab) <- c("n =  50", "n = 100", "n = 200", "n = 300", "n = 400", "n = 500")
tab

df <- data.frame(n = rep(c(50, 100, 200, 300, 400, 500), 4),
                 MSE = unname(c(tab[,1],tab[,2],tab[,3],tab[,4])),
                 type = c( rep("Fixed Thinning ", nrow(tab)), 
                           rep("Random Thinning ", nrow(tab)),
                           rep("Stein Thinning ", nrow(tab)),
                           rep("Full Sample ", nrow(tab))))

library(ggplot2)
ggplot() +
  geom_line(data = df, mapping = aes(x = n, y = MSE, col = type), lwd = 1) +
  labs(x = "Sample Size", y = "MSE") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14))


