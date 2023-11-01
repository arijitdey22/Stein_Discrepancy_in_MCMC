#rm(list = ls())

#===============================================================================
##functions for stein thinning

library(Rcpp)
sourceCpp("00_Stein_Thin_Indices.cpp")

#===============================================================================
##the main function

thin.comp <- function(rep = 1e2, n.all = c(50, 100, 200, 500, 1000)){
  
  #true values
  d = 1
  mu.true <- rep(0,d)
  sigma.true <- diag(d)
  
  #storage for the values
  store.arr <- array(0, dim = c(length(n.all), 4, rep))
  
  #loops start here
  for (i.n in 1:length(n.all))
  {
    for (i.rep in 1:rep)
    {
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
      ind.3 <- stein_thin_indices_cpp(samp, m)
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
source("00_Draw_MVN_MCMC.R")

out <- thin.comp(rep = 25, n.all = c(100, 250, 500, 1000, 2500, 5000, 10000))
save(out, file = "05_Thin_comp_Biased_MCMC.Rdata")

#----------------------
#rm(list = ls())

load("05_Thin_comp_Biased_MCMC.Rdata")
tab <- as.table(round(apply(out, c(1, 2), mean),5))
colnames(tab) <- c("Fixed_Thinning", "Random_Thinning", "Stein_Thinning", "Full Sample")
rownames(tab) <- c("n =  100", "n = 250", "n = 500", "n = 1000", "n = 2500", "n = 5000", "n = 10000")
tab

df <- data.frame(n = rep(c(100, 250, 500, 1000, 2500, 5000, 10000), 4),
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


