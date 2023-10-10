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
  ret <- val^b
  return(ret)
}

#-------------------------------------------------------
#defining the Stein Kernel
k_j0 <- function(x, y, j, b.x, b.y, c = 1){
  
  val <-  c^2 + norm(x-y, type = "2")^2 
  
  term1 <- b.x[j] * b.y[j] * k_xy(x, y, val)
  term2 <- b.x[j] * dy.k_xy(x, y, j, val)
  term3 <- b.y[j] * dx.k_xy(x, y, j, val)
  term4 <- dx.dy.k_xy(x, y, j, val)
  
  return(term1 + term2 + term3 + term4)
  
}

#-------------------------------------------------------
#calculating the pairwise-KSD
Kp.pair <- function(x,y,n)
{
  
  d <- length(x)
  w.sq <- numeric()
  
  b.x <- -x
  b.y <- -y
  
  for ( j in 1:d){
    w.sq[j] <- 1/n^2 * k_j0(x, y, j, b.x, b.y)
  }
  
  ret <- sqrt(sum(w.sq))
  
  return(ret)
}

#-------------------------------------------------------
#choosing the indices

library(mvtnorm)
d <- 5
n <- 100
samp <- rmvnorm(n, mean = rep(0,d), sigma = diag(d))

#calculating KSD of the form K_p(X_i, x_i)
KSD.n <- numeric()
for (i.n in 1:n)
{ 
  KSD.n[i.n] <- Kp.pair(samp[i.n,],samp[i.n,], n)
}

#We are intended to choose the indices
ind.chosen <- numeric()

ind.chosen[1] <- which.min(KSD.n)    #first value is based on KSD.n only

m <- ceiling(log(n)^2)               #setting the value of m              

#loop for j = 2:m
for (j in 2:m){
  
  KSD.trail.n <- numeric()           #storing the required quantity for each n
  
  #loop for n
  for (i in 1:n){          
    f.term <- KSD.n[i] / 2           #first term
    
    s.term <- 0                      #second term
    for (j.d in 1:(j-1))             #loop for j'
    {            
      s.term <- s.term + Kp.pair( samp[ind.chosen[j.d],], samp[i,], n)
    }
    
    KSD.trail.n[i] <- f.term + s.term
  }
  
  ind.chosen[j] <- which.min(KSD.trail.n)
  
  #------------
  #KSD can be computed as we go along
  
  
}














