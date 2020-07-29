



adj.mat.func <- function(y,t){

   N <- length(y);  adj.mat <- matrix(0,nrow = N, ncol = N)

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
   return(adj.mat)        }


wt.mat.func <- function(y, t){  N <- length(y);    wt.mat <- matrix(, nrow = N, ncol = N)
   for(i in 1:N){     for(j in 1:N){ wt.mat[i,j] <- atan( (y[i] - y[j])/(t[i] - t[j]) ) }}
     diag(wt.mat) <- array(0, N) 
      return(wt.mat)         }
                             


