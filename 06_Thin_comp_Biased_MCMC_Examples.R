#rm(list = ls())

#===============================================================================
##functions for stein thinning

library(Rcpp)
sourceCpp("00_Stein_Thin_Indices.cpp")

#--  --  --  --  --  --  --  --  --  --  -- 
#drawing MCMC samples form Normal Distribution

draw.mvn.mcmc <- function(n, mu = 0, sigma = 1, cand.s = 0.5)
{
  
  acp <- 0
  
  ret <- sample(2:5, 1) + rnorm(1)
  
  for (i in 2:n)
  {
    cand.samp <- rnorm (1, ret[i-1], cand.s^2)
    
    HR <- dnorm(cand.samp, mu, sigma) / dnorm(ret[i-1], mu, sigma)
    
    U <- runif(1)
    if (U < HR)
    {
      ret[i] <- cand.samp
      acp <- acp + 1
    }else
    {
      ret[i] <- ret[i-1]
    }
    
  }
  
  #-------------
  
  return(list(samp = ret, acp.rat = (acp/n)))
}

#===============================================================================
##the main function

thin.comp <- function(mu, sigma, n, rep = 1e2){
  
  #true values
  mu.true <- 0
  sigma.true <- 1
  
  #storage for the values
  store.arr <- array(0, dim = c(3, 4, rep))
  store.samp.final <- list()
  
  #loops start here
  for (i.rep in 1:rep)
  {
    store.samp.cur <- list()
    
    set.seed(100 + i.rep)
    
    m <- ceiling(log(n)^2)
    
    samp <- draw.mvn.mcmc(n, mu = mu, sigma = sigma)[[1]]
    
    #---------------------------------------------
    #random thinning
    ind.1 <- sample(1:n, m, replace = T)
    samp.1 <- samp[ind.1]
    mu.est.1 <- mean(samp.1)
    sigma.est.1 <- sd(samp.1)
    
    store.arr[1, 1, i.rep] <- abs(mu.est.1 - mu.true)
    store.arr[2, 1, i.rep] <- abs(qnorm(0.8, mu.est.1, sigma.est.1) - qnorm(0.8))
    store.arr[3, 1, i.rep] <- abs(qnorm(0.9, mu.est.1, sigma.est.1) - qnorm(0.9))
    
    store.samp.cur[[1]] <- samp.1
    
    #---------------------------------------------
    #fixed thinning
    ind.2 <- seq(1, m * floor(n/m), by = floor(n/m))
    samp.2 <- samp[ind.2]
    mu.est.2 <- mean(samp.2)
    sigma.est.2 <- sd(samp.2)
    
    store.arr[1, 2, i.rep] <- abs(mu.est.2 - mu.true)
    store.arr[2, 2, i.rep] <- abs(qnorm(0.8, mu.est.2, sigma.est.2) - qnorm(0.8))
    store.arr[3, 2, i.rep] <- abs(qnorm(0.9, mu.est.2, sigma.est.2) - qnorm(0.9))
    
    store.samp.cur[[2]] <- samp.2
    
    #---------------------------------------------
    #Stein thinning
    ind.3 <- stein_thin_indices_cpp(as.matrix(samp), m)
    samp.3 <- samp[ind.3] 
    mu.est.3 <- mean(samp.3)
    sigma.est.3 <- sd(samp.3)
  
    store.arr[1, 3, i.rep] <- abs(mu.est.3 - mu.true)
    store.arr[2, 3, i.rep] <- abs(qnorm(0.8, mu.est.3, sigma.est.3) - qnorm(0.8))
    store.arr[3, 3, i.rep] <- abs(qnorm(0.9, mu.est.3, sigma.est.3) - qnorm(0.9))
    
    store.samp.cur[[3]] <- samp.3
    
    #---------------------------------------------
    #No thinning
    mu.est.4 <- mean(samp)
    sigma.est.4 <- sd(samp)
    
    store.arr[1, 4, i.rep] <- abs(mu.est.4 - mu.true)
    store.arr[2, 4, i.rep] <- abs(qnorm(0.8, mu.est.4, sigma.est.4) - qnorm(0.8))
    store.arr[3, 4, i.rep] <- abs(qnorm(0.9, mu.est.4, sigma.est.4) - qnorm(0.9))
    
    store.samp.cur[[4]] <- samp
    
    #---------------------------------------------
    store.samp.final[[i.rep]] <- store.samp.cur
    
    #aesthetics
    print("-------------------------------")
    print(paste0("We are at rep = ", i.rep))
    
  }
  
  return(list(values = store.arr, samples = store.samp.final))
  
}

#===============================================================================
###Running the function

n <- 1000

##(a) mu = 0.5, sigma.sq = 1

out.1 <- thin.comp(mu = 0.5, sigma = 1, n, rep = 25)

tab.1 <- as.table(round(apply(out.1[[1]], c(1, 2), mean),5))
samp.all.1 <- out.1[[2]]

## -- -- -- -- 
##(b) mu = 2, sigma.sq = 1

out.2 <- thin.comp(mu = 2, sigma = 1, n, rep = 25)

tab.2 <- as.table(round(apply(out.2[[1]], c(1, 2), mean),5))
samp.all.2 <- out.2[[2]]

## -- -- -- -- 
##(c) mu = 2, sigma.sq = 2

out.3 <- thin.comp(mu = 2, sigma = sqrt(2), n, rep = 25)

tab.3 <- as.table(round(apply(out.3[[1]], c(1, 2), mean),5))
samp.all.3 <- out.3[[2]]

## -- -- -- -- 
##(d) mu = 0, sigma.sq = 5

out.4 <- thin.comp(mu = 0, sigma = sqrt(5), n, rep = 25)

tab.4 <- as.table(round(apply(out.4[[1]], c(1, 2), mean),5))
samp.all.4 <- out.4[[2]]

# save(tab.1, samp.all.1, tab.2, samp.all.2, tab.3, samp.all.3, tab.4, samp.all.4, 
#      file = "06_Thin_comp_Biased_MCMC_Examples.250.Rdata")
# save(tab.1, samp.all.1, tab.2, samp.all.2, tab.3, samp.all.3, tab.4, samp.all.4,
#      file = "06_Thin_comp_Biased_MCMC_Examples.500.Rdata")
# save(tab.1, samp.all.1, tab.2, samp.all.2, tab.3, samp.all.3, tab.4, samp.all.4,
#      file = "06_Thin_comp_Biased_MCMC_Examples.1000.Rdata")
# save(tab.1, samp.all.1, tab.2, samp.all.2, tab.3, samp.all.3, tab.4, samp.all.4, 
#      file = "06_Thin_comp_Biased_MCMC_Examples.5000.Rdata")

#===============================================================================

#rm(list = ls())

# load("06_Thin_comp_Biased_MCMC_Examples.250.Rdata")
# load("06_Thin_comp_Biased_MCMC_Examples.500.Rdata")
# load("06_Thin_comp_Biased_MCMC_Examples.1000.Rdata")
# load("06_Thin_comp_Biased_MCMC_Examples.5000.Rdata")

library(ggplot2)

grid.x <- seq(-6, 10, 0.01)
target.dens <- dnorm(grid.x)

##(a) mu = 0.5, sigma.sq = 1

colnames(tab.1) <- c("Fixed_Thinning", "Random_Thinning", "Stein_Thinning", "Full Sample")
rownames(tab.1) <- c("MSE(Mean)", "MSE(80th Per)", "MSE(90th Per)")
tab.1

#merging all the samples
mer.samp.1 <- vector("list", length = 4)
for (i in 1:4) {
  mer.samp.1[[i]] <- do.call(c, lapply(samp.all.1, function(sublist) sublist[[i]]))
}

set.seed(100 - 1)
samp.one.1 <- samp.all.1[[ sample(1:length(samp.all.1), size = 1) ]]
df.one.1 <- data.frame(samp = c(samp.one.1[[1]], samp.one.1[[2]],
                                samp.one.1[[3]], samp.one.1[[4]]),
                       type = c( rep("Fixed-thinning ", length(samp.one.1[[1]]) ), 
                                 rep("Random-thinning ", length(samp.one.1[[2]]) ),
                                 rep("Stein-thinning ", length(samp.one.1[[3]]) ),
                                 rep("Full Sample ", length(samp.one.1[[4]]) )   )     )

df.all.1 <- data.frame(samp = c(mer.samp.1[[1]], mer.samp.1[[2]],
                                mer.samp.1[[3]], mer.samp.1[[4]]),
                       type = c( rep("Fixed-thinning ", length(mer.samp.1[[1]]) ), 
                                 rep("Random-thinning ", length(mer.samp.1[[2]]) ),
                                 rep("Stein-thinning ", length(mer.samp.1[[3]]) ),
                                 rep("Full Sample ", length(mer.samp.1[[4]]) )   )     )

p.one.1 <- ggplot() +
  geom_density(data = df.one.1, mapping =  aes(x = samp, col = type), lwd = 0.8) + 
  geom_line(aes(x = grid.x, y = target.dens), lwd = 0.8, alpha = 0.65) +
  labs(x = "", y = "Density") + xlim(-3.5,5) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 15))

p.all.1 <- ggplot() +
  geom_density(data = df.all.1, mapping =  aes(x = samp, col = type), lwd = 0.8) + 
  geom_line(aes(x = grid.x, y = target.dens), lwd = 0.8, alpha = 0.65) +
  labs(x = "", y = "Density") + xlim(-3.5,8) + #ylim(0,0.66) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 15))

#---------------------------------------------------------------------------------

##(b) mu = 2, sigma.sq = 1

colnames(tab.2) <- c("Fixed_Thinning", "Random_Thinning", "Stein_Thinning", "Full Sample")
rownames(tab.2) <- c("MSE(Mean)", "MSE(80th Per)", "MSE(90th Per)")
tab.2

#merging all the samples
mer.samp.2 <- vector("list", length = 4)
for (i in 1:4) {
  mer.samp.2[[i]] <- do.call(c, lapply(samp.all.2, function(sublist) sublist[[i]]))
}

set.seed(100 + 2)
samp.one.2 <- samp.all.2[[ sample(1:length(samp.all.2), size = 1) ]]
df.one.2 <- data.frame(samp = c(samp.one.2[[1]], samp.one.2[[2]],
                                samp.one.2[[3]], samp.one.2[[4]]),
                       type = c( rep("Fixed-thinning ", length(samp.one.2[[1]]) ), 
                                 rep("Random-thinning ", length(samp.one.2[[2]]) ),
                                 rep("Stein-thinning ", length(samp.one.2[[3]]) ),
                                 rep("Full Sample ", length(samp.one.2[[4]]) )   )     )

df.all.2 <- data.frame(samp = c(mer.samp.2[[1]], mer.samp.2[[2]],
                                mer.samp.2[[3]], mer.samp.2[[4]]),
                       type = c( rep("Fixed-thinning ", length(mer.samp.2[[1]]) ), 
                                 rep("Random-thinning ", length(mer.samp.2[[2]]) ),
                                 rep("Stein-thinning ", length(mer.samp.2[[3]]) ),
                                 rep("Full Sample ", length(mer.samp.2[[4]]) )   )     )

p.one.2 <- ggplot() +
  geom_density(data = df.one.2, mapping =  aes(x = samp, col = type), lwd = 0.8) + 
  geom_line(aes(x = grid.x, y = target.dens), lwd = 0.8, alpha = 0.65) +
  labs(x = "", y = "Density") + xlim(-3.5,7.5) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 15))

p.all.2 <- ggplot() +
  geom_density(data = df.all.2, mapping =  aes(x = samp, col = type), lwd = 0.8) + 
  geom_line(aes(x = grid.x, y = target.dens), lwd = 0.8, alpha = 0.65) +
  labs(x = "", y = "Density") + xlim(-3.5,8) + #ylim(0,0.66) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 15))

#---------------------------------------------------------------------------------

##(c) mu = 2, sigma.sq = 2

colnames(tab.3) <- c("Fixed_Thinning", "Random_Thinning", "Stein_Thinning", "Full Sample")
rownames(tab.3) <- c("MSE(Mean)", "MSE(80th Per)", "MSE(90th Per)")
tab.3

#merging all the samples
mer.samp.3 <- vector("list", length = 4)
for (i in 1:4) {
  mer.samp.3[[i]] <- do.call(c, lapply(samp.all.3, function(sublist) sublist[[i]]))
}

set.seed(100 - 3)
samp.one.3 <- samp.all.3[[ sample(1:length(samp.all.3), size = 1) ]]
df.one.3 <- data.frame(samp = c(samp.one.3[[1]], samp.one.3[[2]],
                                samp.one.3[[3]], samp.one.3[[4]]),
                       type = c( rep("Fixed-thinning ", length(samp.one.3[[1]]) ), 
                                 rep("Random-thinning ", length(samp.one.3[[2]]) ),
                                 rep("Stein-thinning ", length(samp.one.3[[3]]) ),
                                 rep("Full Sample ", length(samp.one.3[[4]]) )   )     )

df.all.3 <- data.frame(samp = c(mer.samp.3[[1]], mer.samp.3[[2]],
                                mer.samp.3[[3]], mer.samp.3[[4]]),
                       type = c( rep("Fixed-thinning ", length(mer.samp.3[[1]]) ), 
                                 rep("Random-thinning ", length(mer.samp.3[[2]]) ),
                                 rep("Stein-thinning ", length(mer.samp.3[[3]]) ),
                                 rep("Full Sample ", length(mer.samp.3[[4]]) )   )     )

p.one.3 <- ggplot() +
  geom_density(data = df.one.3, mapping =  aes(x = samp, col = type), lwd = 0.8) + 
  geom_line(aes(x = grid.x, y = target.dens), lwd = 0.8, alpha = 0.65) +
  labs(x = "", y = "Density") + xlim(-3.5,6.5) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 15))

p.all.3 <- ggplot() +
  geom_density(data = df.all.3, mapping =  aes(x = samp, col = type), lwd = 0.8) + 
  geom_line(aes(x = grid.x, y = target.dens), lwd = 0.8, alpha = 0.65) +
  labs(x = "", y = "Density") + xlim(-3.5,8) + #ylim(0,0.66) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 15))

#---------------------------------------------------------------------------------

##(c) mu = 2, sigma.sq = 2

colnames(tab.4) <- c("Fixed_Thinning", "Random_Thinning", "Stein_Thinning", "Full Sample")
rownames(tab.4) <- c("MSE(Mean)", "MSE(80th Per)", "MSE(90th Per)")
tab.4

#merging all the samples
mer.samp.4 <- vector("list", length = 4)
for (i in 1:4) {
  mer.samp.4[[i]] <- do.call(c, lapply(samp.all.4, function(sublist) sublist[[i]]))
}

set.seed(100 + 4)
samp.one.4 <- samp.all.4[[ sample(1:length(samp.all.4), size = 1) ]]
df.one.4 <- data.frame(samp = c(samp.one.4[[1]], samp.one.4[[2]],
                                samp.one.4[[3]], samp.one.4[[4]]),
                       type = c( rep("Fixed-thinning ", length(samp.one.4[[1]]) ), 
                                 rep("Random-thinning ", length(samp.one.4[[2]]) ),
                                 rep("Stein-thinning ", length(samp.one.4[[3]]) ),
                                 rep("Full Sample ", length(samp.one.4[[4]]) )   )     )

df.all.4 <- data.frame(samp = c(mer.samp.4[[1]], mer.samp.4[[2]],
                                mer.samp.4[[3]], mer.samp.4[[4]]),
                       type = c( rep("Fixed-thinning ", length(mer.samp.4[[1]]) ), 
                                 rep("Random-thinning ", length(mer.samp.4[[2]]) ),
                                 rep("Stein-thinning ", length(mer.samp.4[[3]]) ),
                                 rep("Full Sample ", length(mer.samp.4[[4]]) )   )     )

p.one.4 <- ggplot() +
  geom_density(data = df.one.4, mapping =  aes(x = samp, col = type), lwd = 0.8) + 
  geom_line(aes(x = grid.x, y = target.dens), lwd = 0.8, alpha = 0.65) +
  labs(x = "", y = "Density") + xlim(-3.5,7) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 15))

p.all.4 <- ggplot() +
  geom_density(data = df.all.4, mapping =  aes(x = samp, col = type), lwd = 0.8) + 
  geom_line(aes(x = grid.x, y = target.dens), lwd = 0.8, alpha = 0.65) +
  labs(x = "", y = "Density") + xlim(-6,10) + #ylim(0,0.66) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 15))


#---------------------------------------------------------------------------------
#all together

library(patchwork)
library(gridExtra)
library(grid)

combined_plot.one <- p.one.1 + p.one.2 + p.one.3 + p.one.4 + plot_layout(ncol = 2) + 
  plot_annotation(theme = theme(legend.position = "bottom")) +
  plot_layout(guides = "collect") & xlab(NULL) & ylab(NULL)

#combined_plot.one

final.plot.one <- wrap_elements(grid::textGrob("Density", rot = 90, x = 0.3, y = 0.47, 
                                  gp = gpar(fontsize = 20))) + combined_plot.one  +
  plot_layout(ncol = 2, widths = c(0.1,2))
final.plot.one

#-- -- -- -- -- 

combined_plot.all <- p.all.1 + p.all.2 + p.all.3 + p.all.4 + plot_layout(ncol = 2) + 
  plot_annotation(theme = theme(legend.position = "bottom")) +
  plot_layout(guides = "collect") & xlab(NULL) & ylab(NULL)

#combined_plot.all

final.plot.all <- wrap_elements(grid::textGrob("Density", rot = 90, x = 0.3, y = 0.47, 
                                               gp = gpar(fontsize = 20))) + combined_plot.all  +
  plot_layout(ncol = 2, widths = c(0.1,2))
#final.plot.all


