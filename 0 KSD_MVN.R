###Defining the function

rm(list = ls())

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
#-------------------------------------------------------
#the final function (two sample)

###inputs
#Samp1: Sample in form of matrix with rows corresponding to different samples
#samp2: A second set of sample with the exact dimension as samp1

# KSD <- function(samp1, samp2)
# {
#   
#   d <- ncol(samp1)
#   n <- nrow(samp1)
#   w <- numeric()
#   store <- numeric()
#   
#   for ( j in 1:d){
#     
#     sum = 0
#     
#     for (i in 1:n){
#       for (i.d in 1:n){
#         
#         xi <- as.vector(samp1[i,])
#         xi.d <- as.vector(samp2[i.d,])
#         
#         Q.n.xi <- Q.n(xi)
#         Q.n.xi.d <- Q.n(xi.d) 
#         k_j0.val <- k_j0(xi, xi.d, j)
#         
#         sum = sum + Q.n.xi * Q.n.xi.d * k_j0.val
#         
#       }
#     }
#     
#     store <- c(store, sum)
#     w[j] <- sqrt(sum)
#     
#     #-------
#     #aesthetics
#     print("----------------------")
#     print(paste0("We are at j = ", j))
#     
#   }
#   
#   ret <- norm(w, type = "2")
#   
#   return(list(ret, store))
# }

#-------------------------------------------------------
#the final function (one sample)

###inputs
#Samp: Sample in form of matrix with rows corresponding to different samples

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
        
        b.x <- b.func(x)
        b.y <- b.func(y)

        sum = sum + 1/n^2 * k_j0(x, y, j, b.x, b.y)
        
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


####################################################################
##Application on Standard Multivariate Normal distribution


#score function (This function is problem specific)
b.func <- function(vec){
  return(-vec)
}

library(mvtnorm)
d <- 5
#samp1 <- rmvnorm(10, mean = rep(0,d), sigma = diag(d))
#samp2 <- rmvnorm(10, mean = rep(0,d), sigma = diag(d))
samp <- rmvnorm(100, mean = rep(0,d), sigma = diag(d))

#KSD(samp1 = samp1, samp2 = samp2)
KSD(samp = samp)


#===================================================================
####################################################################
#===================================================================
##Plot of KSD vs n through parallel computing

rm(list = ls())

#--------------
#defining the function compatible to parallel computing

library(foreach)
library(doParallel)

num_cores <- 6

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

n <- c(10, 25, 50, 100, 200, 500, 1000, 2000, 5000)
KSDs <- numeric()
time.elaps <- numeric()

library(mvtnorm)

for (i in 1:length(n)){
  
  d <- 5
  samp <- rmvnorm(n[i], mean = rep(0,d), sigma = diag(d))
  
  tick <- proc.time()[3]
  
  KSDs[i] <-KSD.par(samp = samp)
  
  tock <- proc.time()[3]
  
  time.elaps[i] <- tock - tick
  
  print(paste0("n = ", n[i]))
}

stopCluster(cl)

save(n, KSDs, time.elaps, file = "0 KSD_MVN_Plot.Rdata")

#--------------

rm(list = ls())
load("0 KSD_MVN_Plot.Rdata")

n
KSDs
time.elaps     #in seconds
#time.elaps.scaled <- time.elaps
time.elaps.scaled <- time.elaps * max(KSDs) / max(time.elaps)   #scaling of the time component

library(ggplot2)

df <- data.frame(c(n,n),
                 c(KSDs,time.elaps.scaled),
                 c(rep("KSD", length(KSDs)), rep("Time", length(time.elaps.scaled)) ) )
colnames(df) <- c("n", "Values", "Legend")

ggplot() +
  geom_line(data = df, mapping = aes(x = n, y = Values, color = Legend), lwd = 1) +
  labs(x = "Sample size", y = "KSD",
       title = "'Sample size vs KSD' and 'Sample size vs Time' Plot",
       caption = "Note: For optimal visualization, time component has been scaled
       according as the KSD component") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        plot.caption = element_text(size = 12))
