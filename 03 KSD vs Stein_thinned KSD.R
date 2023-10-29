#rm(list = ls())

#-------------------------------------------------------------------------------
#defining the Stein Kernel
k_j0 <- function(x, y, j, b.x, b.y, b = -0.5, c = 1){
  
  val <-  c^2 + norm(x-y, type = "2")^2 
  
  term1 <- b.x[j] * b.y[j] * (  val^b  )
  term2 <- b.x[j] * (  b * val^(b-1) * (-2) * (x[j] - y[j])  )
  term3 <- b.y[j] * (  b * val^(b-1) *   2  * (x[j] - y[j])  )
  term4 <- (  - 2 * b * (  (b-1) * val^(b-2) * 2 * (x[j] - y[j])^2  +  val^(b-1)  )  )
  
  return(term1 + term2 + term3 + term4)
  
}

#===============================================================================
#calculating the pairwise-KSD
Kp.xy <- function(x,y,n)
{
  
  d <- length(x)
  w.sq <- numeric()
  
  b.x <- -x
  b.y <- -y
  
  for ( j in 1:d){
    w.sq[j] <-  k_j0(x, y, j, b.x, b.y)
  }
  
  ret <- sum(w.sq)
  
  return(ret)
}

#-------------------------------------------------------
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
#function for calculating the KSD

KSD <- function(samp)
{
  
  d <- ncol(samp)
  n <- nrow(samp)
  w <- numeric()
  
  for ( j in 1:d){
    
    sum = 0
    
    for (i.x in 1:n){
      for (i.y in 1:n){
        
        x <- as.vector(samp[i.x,])
        y <- as.vector(samp[i.y,])
        
        b.x <- -x
        b.y <- -y
        
        sum = sum + 1/n^2 * k_j0(x, y, j, b.x, b.y)
        #sum = c(sum, 1/n^2 * k_j0(x, y, j, b.x, b.y))
        
      }
    }
    
    w[j] <- sqrt(sum)
    
    #-------
    #aesthetics
    # print("----------------------")
    # print(paste0("We are at j = ", j))
    
  }
  
  ret <- norm(w, type = "2")
  
  return(ret)
}

#===============================================================================
#calculating the KSD vs thinned KSD

source("00 Draw_MVN_MCMC.R")

n.all <- c(100, 150, 200, 250, 300, 400, 500, 1000, 2500, 5000)
KSD.all <- numeric()
time.elaps.KSD <- numeric()

KSD.thinned <- numeric()
time.elaps.KSD.thinned <- numeric()
time.elaps.KSD.thinned.total <- numeric()

time.elaps.KSD.thinned1 <- numeric()

d <- 2

for (i in 1:length(n.all)){
  
  set.seed(100 + i)
  
  n <- n.all[i]
  samp <- draw.mvn.mcmc(n, d = d, mu = rep(0,d), sigma = diag(d))[[1]]
  
  # -- -- --
  
  tick <- proc.time()[3]
  KSD.all[i] <- KSD(samp = samp)
  tock <- proc.time()[3]
  time.elaps.KSD[i] <- tock - tick
  
  # -- -- -- 
  
  tick <- proc.time()[3]
  m <- ceiling(log(n)^2)
  ind <- stein.thin.indices(samp, m)
  samp.thin <- samp[ind, ]
  
  tick1 <- proc.time()[3]
  KSD.thinned[i] <- KSD(samp = samp.thin)
  tock <- proc.time()[3]
  
  # -- -- -- 
  
  time.elaps.KSD.thinned[i] <- tock - tick1
  time.elaps.KSD.thinned1[i] <- tock - tick1
  time.elaps.KSD.thinned.total[i] <- tock - tick
  
  # -- -- -- 
  print(paste0("We are at n = ", n.all[i]))
  
}

save(n.all, KSD.all, time.elaps.KSD, KSD.thinned, time.elaps.KSD.thinned, time.elaps.KSD.thinned.total,
     file = "03 KSD vs Stein_thinned KSD.Rdata")


#===============================================================================

#rm(list = ls())
load("03 KSD vs Stein_thinned KSD.Rdata")

df <- data.frame(n = c(n.all,n.all), 
                 KSD = c(KSD.all, KSD.thinned),
                 type = c(rep( "Full Sample", length(KSD.all) ), 
                          rep("Thinned Sample", length(KSD.thinned) )))

library(ggplot2)
ggplot() +
  geom_line(data = df, mapping = aes(x = n, y = log(KSD), col = type), lwd = 1) +
  labs(x = "Sample size", y = "KSD") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 16)) + xlim(c(100,2500))

#--------------------------------------------
#plot for time

df.1 <- data.frame(n = rep(n.all,2),
                 Event = c(rep("event1", length(n.all)), 
                           rep("event2", length(n.all))),
                 Time = c(log(time.elaps.KSD),
                          log(time.elaps.KSD.thinned)))

ggplot() +
  geom_line(data = df.1, mapping = aes(x = n, y = Time, col = Event), lwd = 1) +
  labs(x = "Sample Size", y = "log(Time)", color = "") + xlim(c(100,2500)) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        legend.position = "none")

#-- -- -- -- 
#all three together

df.2 <- data.frame(n = rep(n.all,3),
                 Event = c(rep("event1", length(n.all)), 
                           rep("event2", length(n.all)),
                           rep("event3", length(n.all))),
                 Time = c(time.elaps.KSD,
                          time.elaps.KSD.thinned,
                          time.elaps.KSD.thinned.total))

#-- -- -- -- -- --

ggplot() +
  geom_line(aes(x = n.all, y = time.elaps.KSD), color = "#F8766D", size = 1) + 
  geom_line(aes(x = n.all, y = time.elaps.KSD.thinned), color = "#00BFC4", lty = "twodash", size = 1) +
  geom_line(aes(x = n.all, y = time.elaps.KSD.thinned.total), color = "#00BFC4", size = 1) +
  labs(x = "Sample Size", y = "Time") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14))







