#This is the function to calculate balance
#X is the relative abundance matrix
#z is the group indicator vector; 
#  allowing NA values in z corresponds to undetermined group status

balance <- function(X,z){
  if (is.null(dim(X))){X <- matrix(X,nrow=1)}
  z[is.na(z)] <- 0
  bplus <- z == 1 
  mplus <- sum(bplus)
  bminus <- z == -1
  mminus <- sum(bminus)
  
  if (mplus == 0 && mminus == 0) {B<- rep(0,dim(X)[1])}
  if (mplus == 0 && mminus != 0){
    B <- apply(X,1,
               function(x){
                 x <- log(x) # take log
                 - sum(x[bminus]) / sqrt(mminus) 
               })
  }
  if (mminus == 0 && mplus != 0){
    B <- apply(X,1,
               function(x){
                 x <- log(x) # take log
                 sum(x[bplus]) / sqrt(mplus) 
               })
  } 
  if (mminus !=0 && mplus !=0){
    B <- apply(X,1,
               function(x){
                 x <- log(x) # take log
                 (sum(x[bplus]) / mplus - sum(x[bminus]) / mminus) *
                   sqrt(1 / (1/mplus + 1/mminus))
               })
  }
  
  
  return(B)
}


alr <- function(P){
  q <- dim(P)[2]
  t(apply(P,1,function(x){log(x[-q]/x[q])}))
}
alr_inv <- function(alrP){
  #input is a row vector in matrix format
  n <- dim(alrP)[1]
  alrP <- exp(cbind(alrP,0))
  return(matrix(mapply('/',alrP,rowSums(alrP)),nrow=n,byrow=F))
}