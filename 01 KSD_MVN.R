#rm(list = ls())

##Necessary functions
#-------------------------------------------------------
dx.dy.k_xy <- function(x, y, j, val, b = -0.5, c = 1){
  ret <- - 2 * b * (  (b-1) * val^(b-2) * 2 * (x[j] - y[j])^2  +
                      val^(b-1)  )
  return(ret)
}

#-------------------------------------------------------
dx.k_xy <- function(x, y, j, val, b = -0.5, c = 1){
  ret <- b * val^(b-1) * 2 * (x[j] - y[j])
  return(ret)
}

#-------------------------------------------------------
#
dy.k_xy <- function(x, y, j, val, b = -0.5, c = 1){
  ret <- b * val^(b-1) * (-2) * (x[j] - y[j])
  return(ret)
}

#-------------------------------------------------------
#The KSD with IMQ base kernel
k_xy <- function(x, y, val, b = -0.5, c = 1){
  if (length(x) == length(y)){
    ret <- val^b
  }
  else{
    ret <- Inf
  }
  return(ret)
}

#-------------------------------------------------------
#defining the Stein Kernel
k_j0 <- function(x, y, j, b.x, b.y, c = 1){
  
  val <-  c^2 + norm(x-y, type = "2")^2 
  
  term1 <- b.x[j] * b.y[j] * k_xy(x,y, val)
  term2 <- b.x[j] * dy.k_xy(x, y, j, val)
  term3 <- b.y[j] * dx.k_xy(x, y, j, val)
  term4 <- dx.dy.k_xy(x, y, j, val)
  
  return(term1 + term2 + term3 + term4)
  
}

#-------------------------------------------------------
#the final function

###input: Samp: Sample in form of matrix with rows corresponding to different samples

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
    print("----------------------")
    print(paste0("We are at j = ", j))
    
  }
  
  ret <- norm(w, type = "2")
  
  return(ret)
}


#===================================================================
####################################################################
##calculation of KSD vs n through parallel computing


#--------------
#defining the function compatible to parallel computing

library(foreach)
library(doParallel)

num_cores <- detectCores() - 2

cl <- makeCluster(num_cores)
registerDoParallel(cl)

KSD.par <- function(samp)
{
  
  parallel_function <- function(j){
    
    sum = 0
    
    for (i.x in 1:n){
      for (i.y in 1:n){
        
        x <- as.vector(samp[i.x,])
        y <- as.vector(samp[i.y,])
        
        b <- -0.5
        c <- 1 
        val <-  c^2 + norm(x-y, type = "2")^2
        
        #------------------
        #b function here
        b.x <- -x
        b.y <- -y
        
        #------------------
        #k_xy function here
        k_xy <- val^b
        
        #------------------
        #dx.dy.k_xy function here
        dx.dy.k_xy <- - 2 * b * (  (b-1) * val^(b-2) * 2 * (x[j] - y[j])^2 + val^(b-1)  )
        
        #------------------
        #dx.k_xy function here
        dx.k_xy <- b * val^(b-1) * 2 * (x[j] - y[j])
        
        #------------------
        #dy.k_xy function here
        dy.k_xy <- b * val^(b-1) * (-2) * (x[j] - y[j])
        
        #------------------
        #k_j0 function here
        term1 <- b.x[j] * b.y[j] * k_xy
        term2 <- b.x[j] * dy.k_xy
        term3 <- b.y[j] * dx.k_xy
        term4 <- dx.dy.k_xy
        
        term = term1 + term2 + term3 + term4
        #------------------
        
        sum = sum + 1/n^2 * term
        
      }
    }
    
    return(sqrt(sum))
    
  }
  
  #-- -- -- -- 
  
  d <- ncol(samp)
  n <- nrow(samp)
  w <- numeric()
  
  list.w <- foreach(j = 1:d) %dopar% {
    parallel_function(j)
  }
  
  for (ii in 1:length(list.w))
  {
    w[ii] <- list.w[[ii]]
  }
  
  ret <- norm(w, type = "2")
  
  return(ret)
}

#--------------
#running the function

n <- c(25, 50, 100, 200, 500, 1000, 2000)
d <- 2
KSDs <- numeric()
time.elaps <- numeric()

source("00 Draw_MVN_MCMC.R")

for (i in 1:length(n)){
  
  set.seed(100+i)
  
  samp <- draw.mvn.mcmc(n[i], d = d, mu = rep(0,d), sigma = diag(d))[[1]]
  
  tick <- proc.time()[3]
  
  KSDs[i] <-KSD.par(samp = samp)
  
  tock <- proc.time()[3]
  
  time.elaps[i] <- tock - tick
  
  print(paste0("n = ", n[i]))
}

stopCluster(cl)

save(n, KSDs, time.elaps, file = "01 KSD_MVN_Plot.Rdata")

#--------------

#rm(list = ls())
load("01 KSD_MVN_Plot.Rdata")

(n = n)
(KSD = KSDs)
time.elaps     #in seconds
time.elaps.scaled <- time.elaps * max(KSDs) / max(time.elaps)   #scaling of the time component
Time <- time.elaps * max(KSDs) / max(time.elaps)   #scaling of the time component

library(ggplot2)

ggplot()+
  geom_line(aes(n,KSD, col = 'KSD'), lwd = 1)+
  geom_point(aes(n,KSD, col = 'KSD'), cex = 1.5, pch = 15)+
  geom_line(aes(n, Time, col = 'Time'), lwd = 1)+
  geom_point(aes(n,Time, col = 'Time'), cex = 1.5, pch = 15)+
  scale_y_continuous(name = "KSD", sec.axis = sec_axis(~./(max(KSDs) / max(time.elaps)), name="Time")) +
  labs(x = "Sample size", y = "KSD") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.title.y.left=element_text(color="#F8766D"),
        axis.text.y.left=element_text(color="#F8766D"),
        axis.title.y.right=element_text(color="#00BFC4"),
        axis.text.y.right=element_text(color="#00BFC4"),
        legend.position = "bottom")
