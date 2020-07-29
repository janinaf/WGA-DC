





wt.vga.mat.func <- function(y,t){

   N <- length(y);  adj.mat <- matrix(0,nrow = N, ncol = N)
                     wt.mat <- matrix(, nrow = N, ncol = N)

   for(i in 1:N){     for(j in 1:N){
        if(abs(j-i) <= 1){adj.mat[i,j] <- 1}
        if(abs(j-i) > 1){
           if(j > i){
             y.temp <- y[(i+1):(j-1)]
             t.temp <- t[(i+1):(j-1)]
             n <- length(t.temp)
             ind <- array(0,n)
              for(k in 1:n){
       xx <- y.temp[k] - (y[j] + ((y[i] - y[j])*((t[j] - t.temp[k])/(t[j] - t[i]))))
              if(xx < 0){ind[k] <- 1}
                           }
             sum.ind <- sum(ind)
             if(sum.ind == n){ adj.mat[i,j] <- 1}}}}}

   adj.mat <- adj.mat + t(adj.mat)
   adj.mat[adj.mat == 2] <- 1

   for(i in 1:N){     for(j in 1:N){ wt.mat[i,j] <- atan( (y[i] - y[j])/(t[i] - t[j]) ) }}
     diag(wt.mat) <- array(0, N)

         wvga.mat <- adj.mat*wt.mat
         return(wvga.mat)         }






sim.normal.dat  <- read.table(".../sim.data.d2a.normal.csv", header = TRUE, sep = ",")
sim.cauchy.dat  <- read.table(".../sim.data.d2a.cauchy.csv", header = TRUE, sep = ",")



source(".../various.functions.r")

library(MASS)
library(psych)
library(ccgarch)
library(rmgarch)
library(fMultivar)

T <- 600 # Number of random samples
#set.seed(123)

d2a.normal.dat  <- read.table(".../normal.d2a.corr.csv", header = FALSE, sep = ",")
d2a.cauchy.dat  <- read.table(".../cauchy.d2a.corr.csv", header = FALSE, sep = ",")


dnm <- dcc.normal.mat <- d2a.normal.dat[,(2*c(0:999)) + 1]
dcm <- dcc.cauchy.mat <- d2a.cauchy.dat[,(2*c(0:999)) + 1]


dcc.mat <- dcm ## change this between cauchy and normal
sim.dat <- sim.cauchy.dat ## change this between cauchy and normal

sim.num <- dim(dnm)[2]


mean.abs.sw.cor.vec <- mean.abs.vga.cor.vec <- mean.abs.dcc.cor.vec <- array(, sim.num)
max.abs.sw.cor.vec <- max.abs.vga.cor.vec <- max.abs.dcc.cor.vec <- array(, sim.num)
sw.cor.mse.vec <- vga.cor.mse.vec <- dcc.cor.mse.vec <- array(, sim.num)


rho.vec <- array(,T)
DDelta <- 1024/(2^3)

for(t in 1:T){ rho.vec[t] <- rho <- sin(t/DDelta) }


for(kk in 1:sim.num){


ts.1 <- sim.dat[,(2*kk)-1]
ts.2 <- sim.dat[,(2*kk)]


dcc.corr <- dcc.mat[,kk]

## new method

ws <- 15 # window.size

sw.corr <- array(,length(ws:T)) #sliding window correlation (initialization)
sim.vec <- array(,length(ws:T)) #similarity vector based on VGA (initialization)

for(i in ws:T){ sw.corr[i] <- cor.test(ts.1[(i-ws+1):i], ts.2[(i-ws+1):i])$est 
                W.1 <- wt.mat.func(ts.1[(i-ws+1):i], c(1:ws))
                W.2 <- wt.mat.func(ts.2[(i-ws+1):i], c(1:ws))
                w.1 <- W.1[upper.tri(W.1)]; w.2 <- W.2[upper.tri(W.2)]
                sim.vec[i] <- cor.test(w.1, w.2)$est }  

     W.1.all <- wt.mat.func(ts.1, c(1:T))
     W.2.all <- wt.mat.func(ts.2, c(1:T))


sim.all.vec <- array(,T) #similarity vector based on VGA (initialization)
sim.all.vec1 <- array(,length(ws:T)) #similarity vector based on VGA (initialization)



for(i in ws:T){ x1 <- apply(W.1.all[(i-ws+1):i,], 2, median); x2 <- apply(W.2.all[(i-ws+1):i,],2,median)
       sim.all.vec1[(i)] <- cor.test(x1,x2)$est  } 



mean.abs.sw.cor.vec[kk]    <- mean(abs(sw.corr), na.rm = T)
mean.abs.vga.cor.vec[kk]   <- mean(abs(sim.all.vec1), na.rm = T)
mean.abs.dcc.cor.vec[kk]  <- mean(abs(dcc.corr), na.rm = T)

max.abs.sw.cor.vec[kk]    <- max(abs(sw.corr), na.rm = T)
max.abs.vga.cor.vec[kk]   <- max(abs(sim.all.vec1), na.rm = T)
max.abs.dcc.cor.vec[kk]  <- max(abs(dcc.corr), na.rm = T)


sw.cor.mse.vec[kk]   <- mean( (sw.corr - rho.vec)^2, na.rm = T)
vga.cor.mse.vec[kk]  <- mean( (sim.all.vec1 - rho.vec)^2, na.rm = T)
dcc.cor.mse.vec[kk] <- mean( (dcc.corr - rho.vec)^2)


print(kk)}




mean(sw.cor.mse.vec)
mean(vga.cor.mse.vec)
mean(dcc.cor.mse.vec)

quantile(sw.cor.mse.vec)
quantile(vga.cor.mse.vec)
quantile(dcc.cor.mse.vec)




indx <- c(rep(1, sim.num), rep(2, sim.num), rep(3, sim.num))
rslts <- c(sw.cor.mse.vec, vga.cor.mse.vec, dcc.cor.mse.vec)


#boxplot(rslts[!is.na(rslts)]~indx[!is.na(rslts)], col = c("green", "magenta", "red"), ylim = c(min(rslts[!is.na(rslts)]),max(rslts[!is.na(rslts)])), xaxt = "n")
boxplot(rslts~indx, col = c("green", "magenta", "red"), ylim = c(min(rslts[!is.na(rslts)]),max(rslts[!is.na(rslts)])), xaxt = "n")
mtext("SW",   side=1, line=1, at=1, cex=1.5, las=1, font=1)
mtext("WGA", side=1, line=1, at=2, cex=1.5, las=1, font=1)
mtext("DCC", side=1, line=1, at=3, cex=1.5, las=1, font=1)





mean(mean.abs.sw.cor.vec)
mean(mean.abs.vga.cor.vec)
mean(mean.abs.dcc.cor.vec)


sd(mean.abs.sw.cor.vec)
sd(mean.abs.vga.cor.vec)
sd(mean.abs.dcc.cor.vec)


mean(max.abs.sw.cor.vec)
mean(max.abs.vga.cor.vec)
mean(max.abs.dcc.cor.vec)


sd(max.abs.sw.cor.vec)
sd(max.abs.vga.cor.vec)
sd(max.abs.dcc.cor.vec)





#plot(sw.corr, type = "l", col = "green", ylab = "dynamic correlation", ylim = c(-1, 1), xlab = "time", cex.axis = 2.5, cex.lab = 2.5, cex.main = 2.5)
plot(sw.corr, type = "l", col = "green", ylab = "", ylim = c(-1, 1), xlab = "", cex.axis = 2.5, cex.lab = 2.5, cex.main = 2.5)
lines(sim.all.vec1, col = "magenta")
lines(dcc.corr, col = "red")
lines(rho.vec, lwd = 2)




