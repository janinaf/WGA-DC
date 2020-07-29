


sim.num <- 10000


source("../various.functions.r")

library(MASS)
library(psych)
#library(ccgarch)
#library(rmgarch)
library(fMultivar)

T <- 1000 # Number of random samples
#set.seed(123)



for(kk in 1:sim.num){


bvn.ts <- matrix(, nrow = T, ncol = 2)
rho.vec <- array(,T)
DDelta <- 1024/(2^4)


for(t in 1:T){

rho.vec[t] <- rho <- 0
mu1 <- 0; s1 <- sqrt(2)
mu2 <- 0; s2 <- sqrt(3)

mu <- c(mu1,mu2) # Mean 
sigma <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2),
           2) # Covariance matrix

#bvn.ts[t,] <- mvrnorm(1, mu = mu, Sigma = sigma ) 
bvn.ts[t,] <- rcauchy2d(1, rho = rho) 
              }
colnames(bvn.ts) <- c("ts.1","ts.2")

ts.1a <- bvn.ts[,1]
ts.2 <- bvn.ts[,2]

thresh <- th <- 50
ts.1a[ts.1a >  th] <-  th
ts.1a[ts.1a < -th] <- -th
ts.2[ts.2 >  th] <-  th
ts.2[ts.2 < -th] <- -th


xx <- array(, T) #vector("numeric",T)
xx[1:3] <- ts.1a[1:3]
#for(t in 4:T){ xx[t] <- 1.3*xx[t-1] + 1.55*xx[t-3] + ts.1a[t] }
for (i in 3:T){ xx[i] <- 0.5*xx[i - 1] + xx[i - 1] - 0.75*xx[i - 2] + ts.1a[i] + 0.3*ts.1a[i - 1]}
#ts.1 <- xx
ts.1 <- ts.1a

ts.dat <- cbind(ts.1, ts.2)

  if(kk == 1){ sim.dat <- ts.dat }
  if(kk > 1){  sim.dat <- cbind(sim.dat, ts.dat) }

 print(kk) }

 
#write.table(sim.dat, file = "/sim.data.d1d.cauchy.csv", sep = ",", row.names = FALSE, na = ".")


