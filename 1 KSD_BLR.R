rm(list = ls())

titanic <- read.csv("https://dvats.github.io/assets/titanic.csv")
Y <- titanic[,1]
X <- as.matrix(titanic[, -1])

head(Y)
head(X)

#---------------------------------------------------
#running the MCMC chain

library(rjags)

mod_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dbern(p[i])
        logit(p[i]) = b[1]*intercept[i] + b[2]*Is_male[i] + b[3]*Age[i] + b[4]*Sib.Sp[i] + b[5]*Par.Ch[i] + b[6]*fare[i]
    }

    for (j in 1:6) {
        b[j] ~ dnorm(0.0, 1.0/100.0)
    }
} "

set.seed(100)
data_jags = list(y = Y, intercept = X[,1], Is_male = X[,2],
                 Age = X[,3], Sib.Sp = X[,4], Par.Ch = X[,5], fare = X[,6])

params = c("b")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=1)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                        variable.names=params,
                        n.iter=5e3)

samp <- as.matrix(mod_sim)
dim(samp)

#---------------------------------------------------
#function needed for computing KSD

library(foreach)
library(doParallel)

num_cores <- 6

cl <- makeCluster(num_cores)
registerDoParallel(cl)

KSD.par <- function(samp, X, Y)
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
        #b function here (this is only the problem specific part)
        beta <- matrix(x, ncol = 1)
        p_i <- 1/(1+exp(-X %*% beta))
        b.x <- -2*x + colSums(X * as.numeric(Y - p_i))
        
        beta <- matrix(y, ncol = 1)
        p_i <- 1/(1+exp(-X %*% beta))
        b.y <- -2*y + colSums(X * as.numeric(Y - p_i))
        
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

#---------------------------------------------------
#running the function

#KSD.par(samp = samp)
KSD.par(samp = samp[1:100, ], X, Y)

#===================================================================
##Plot of KSD vs n through parallel computing

n <- c(10, 25, 50, 100, 200, 500, 1000, 2000)
KSDs <- numeric()
time.elaps <- numeric()

for (i in 1:length(n)){
  
  samp.loop <- samp[1:n[i], ]
  
  tick <- proc.time()[3]
  
  KSDs[i] <-KSD.par(samp = samp.loop, X, Y)
  
  tock <- proc.time()[3]
  
  time.elaps[i] <- tock - tick
  
  print(paste0("n = ", n[i]))
  
  
}

stopCluster(cl)

save(n, KSDs, time.elaps, file = "1 KSD_BLR_Plot.Rdata")

#--------------

rm(list = ls())
load("1 KSD_BLR_Plot.Rdata")

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

