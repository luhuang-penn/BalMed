#Metroplis Hastings Algorithm
#proposal density:
# phi0 probability for changes between 0 and 1 (addition,deletion and switch)
# phi1 probability for changes between 0 and -1
# phi2 probaiblity for changes between 1 and -1
# pa, probability for addition
# pd,probability for deletion
# pw, probability for swtich

#return is z for next iteration; logtarget value; MH ratio; posterior for beta|z
#possible to add simulated annealing


logtarget <- function(z,Y,P,Trt,hyperparam){

  B <- balance(P,z)
  mplus <- sum(z == 1)
  mminus <- sum(z == -1)
  w1 <- hyperparam$w1
  w2 <- hyperparam$w2
  v <- hyperparam$v
  lambda <- hyperparam$lambda
  n <- dim(P)[1]
  q <- dim(P)[2]
  B1 <- cbind(1,Trt,B)
  beta0 <- hyperparam$beta0
  V <- hyperparam$V
  #b0_ <- matrix(c(a0,b0),ncol = 1)
  #V <- diag(c(h,varc))
  V_star_inv <- solve(V)+t(B1)%*% B1
  V_star <- solve(V_star_inv)
  u <- V_star %*% (solve(V)%*%beta0 + t(B1)%*% Y)
  b_star <- c(v*lambda+t(beta0)%*% V_star_inv %*% beta0 + t(Y) %*% Y - t(u) %*% V_star_inv %*% u)
  logtarget <- log(det(V_star_inv))*(-0.5) - 
    (v + n) / 2 * log(b_star/2) +
    mplus * log(w1) + log(w2) *mminus +  log(1- w1 - w2)*(q - mplus - mminus)
  return(list(logtarget=logtarget,post_beta=list(mean=u,sigma = b_star/(lambda + n)*V_star)))
}

proposal_z <- function(zt,phi0,phi1,pa,pd,pw,verbose){
  znew <- zt
  phi2 <- 1 - phi0 - phi1
  mplus <- sum(zt == 1)
  mminus <- sum(zt == -1)
  m0 <- length(zt) - mplus - mminus
  logratio_proposal <- 0
  #type[1] between 0 and 1
  #type[2] between 0 and -1
  #type[3] between 1 and -1
  type <- rmultinom(1,1,c(phi0,phi1,phi2))
  if (m0 == 0 && (mplus == 1 && mminus == 1)){
    type <- c(0,0,1)
  } 
  while (type[1] && m0 == 0 && mplus == 1 && mminus > 1){
    type <- rmultinom(1,1,c(phi0,phi1,phi2))
  }
  while (type[2] && m0 == 0 && mplus > 1 && mminus == 1){
    type <- rmultinom(1,1,c(phi0,phi1,phi2))
  }
  #between 0 and 1
  #movetype == 1 : 0->1
  #movetype == 2 : 1->0
  #movetype == 3 : 0 <==>1 switch
  
  if (type[1]){
    movetype <- which(rmultinom(1,1,c(pa,pd,pw))==1)
    #when there's no 0 in zt
    if ( m0 == 0 && mplus > 1){
      movetype <- 2
    }
    #when there's only one 1 in zt and no zero
    while (movetype == 2  && mplus == 1 && m0 > 0){
      movetype <- which(rmultinom(1,1,c(pa,pd,pw))==1)
    }
    #movetype == 1 : 0->1
    if (movetype == 1){
      idx <- which(zt == 0)
      ad <-  idx[ceiling(runif(1,0,m0))]
      znew[ad] <- 1
      logratio_proposal <- log(pd/(mplus+1)) - 
        log(pa/m0)
      if (verbose) {
        cat(sprintf('  change 0 to 1 at position %s \n', ad))
      }
    }
    #movetype == 2 : 1->0
    if (movetype == 2){
      idx <- which(zt == 1)
      rm <-  idx[ceiling(runif(1,0,mplus))]
      znew[rm] <- 0
      logratio_proposal <- log(pa/(m0+1)) - 
        log(pd/mplus)
      if (verbose) {
        cat(sprintf('  change 1 to 0 at position %s \n', rm))
      }
    }
    #movetype == 3: 1<=>0 swtich
    #logratio_proposal does not change
    if (movetype == 3){
      idx1 <- which(zt == 0)
      idx2 <- which(zt == 1)
      sw1 <- idx1[ceiling(runif(1,0,m0))]
      sw2 <- idx2[ceiling(runif(1,0,mplus))]
      znew[sw1] <- 1
      znew[sw2] <- 0
      if (verbose) {
        cat(sprintf('  switch 0 and 1 bw pos %s and pos %s \n', sw1,sw2))
      }
    }
    
    
  }
  
  
  #between 0 and -1
  #movetype == 1 : 0->-1
  #movetype == 2 : -1->0
  #movetype == 3 : 0 <==>-1 switch
  if (type[2]){
    movetype <- which(rmultinom(1,1,c(pa,pd,pw))==1)
    #when there's no 0 in zt
    if ( m0 == 0 && mminus > 1){
      movetype <- 2
    }
    #when there's only one -1 in zt
    while (movetype == 2  && mminus == 1 && m0 > 0){
      movetype <- which(rmultinom(1,1,c(pa,pd,pw))==1)
    }
    #movetype == 1 : 0->-1
    if (movetype == 1){
      idx <- which(zt == 0)
      ad <-  idx[ceiling(runif(1,0,m0))]
      znew[ad] <- -1
      logratio_proposal <- log(pd/(mminus+1)) - 
        log(pa/m0)
      if (verbose) {
        cat(sprintf('  change 0 to -1 at position %s \n', ad))
      }
    }
    #movetype == 2 : -1->0
    if (movetype == 2){
      idx <- which(zt == -1)
      rm <-  idx[ceiling(runif(1,0,mminus))]
      znew[rm] <- 0
      logratio_proposal <- log(pa/(m0+1)) - 
        log(pd/mminus)
      if (verbose) {
        cat(sprintf('  change -1 to 0 at position %s \n', rm))
      }
    }
    #movetype == 3: -1<=>0 swtich
    #logratio_proposal does not change
    if (movetype == 3){
      idx1 <- which(zt == 0)
      idx2 <- which(zt == -1)
      sw1 <- idx1[ceiling(runif(1,0,m0))]
      sw2 <- idx2[ceiling(runif(1,0,mminus))]
      znew[sw1] <- -1
      znew[sw2] <- 0
      if (verbose) {
        cat(sprintf('  switch -1 and 0 bw pos %s and pos %s \n', sw1,sw2))
      }
    }
    
  }
  #between 1 and -1
  #movetype == 1 : 1 -> -1
  #movetype == 2 : -1 -> 1
  #movetype == 3 : 1 <==> -1 switch
  if (type[3]){
    movetype <- which(rmultinom(1,1,c(pa,pd,pw))==1)
    
    #when there's only one -1 in zt
    while (movetype == 2  && mminus == 1 && mplus > 1){
      movetype <- which(rmultinom(1,1,c(pa,pd,pw))==1)
    }
    #when there's only one 1 in zt
    while (movetype == 1  && mplus == 1 && mminus > 1){
      movetype <- which(rmultinom(1,1,c(pa,pd,pw))==1)
    }
    #when there's only one 1 and only one -1
    if (mplus == 1 && mminus ==1){
      movetype <- 3
    }
    #movetype == 1 : 1 -> -1
    if (movetype == 1){
      idx <- which(zt == 1)
      ad <-  idx[ceiling(runif(1,0,mplus))]
      znew[ad] <- -1
      logratio_proposal <- log(pd/(mminus+1)) - 
        log(pa/mplus)
      if (verbose) {
        cat(sprintf('  change 1 to -1 at pos %s \n', ad))
      }
    }
    #movetype == 2 : -1->1 pd
    if (movetype == 2){
      idx <- which(zt == -1)
      rm <-  idx[ceiling(runif(1,0,mminus))]
      znew[rm] <- 1
      logratio_proposal <- log(pa/(mplus+1)) - 
        log(pd/mminus)
      if (verbose) {
        cat(sprintf('  change -1 to 1 at pos %s \n', rm))
      }
    }
    #movetype == 3: 1<=> -1 swtich
    #logratio_proposal does not change
    if (movetype == 3){
      idx1 <- which(zt == 1)
      idx2 <- which(zt == -1)
      sw1 <- idx1[ceiling(runif(1,0,mplus))]
      sw2 <- idx2[ceiling(runif(1,0,mminus))]
      znew[sw1] <- -1
      znew[sw2] <- 1
      if (verbose) {
        cat(sprintf('  switch 1 and -1 bw pos %s and pos %s \n', sw1,sw2))
      }
    }
    
  }
  return(list(znew=znew,logratio_proposal = logratio_proposal))
}

iteration_z <- function(Y,P,Trt,zt,hyperparam,phi0=1/3,phi1=1/3,pa=1/3,pd=1/3,pw=1/3,verbose=FALSE){
  #print(dim(X))
  out <- proposal_z(zt,phi0,phi1,pa,pd,pw,verbose)
  znew <- out$znew
  logratio_proposal <- out$logratio_proposal
  #Force regression coefficient of balance to >0
  B <- balance(P,znew)
  fit <- lm(Y~Trt+B)
  while (fit$coefficients[3]<0){
    out <- proposal_z(zt,phi0,phi1,pa,pd,pw,verbose)
    znew <- out$znew
    logratio_proposal <- out$logratio_proposal
  #Force regression coefficient of balance to >0
    B <- balance(P,znew)
    fit <- lm(Y~Trt+B)
  }
  
  u <- runif(1,0,1)
  newltarget <- logtarget(znew,Y,P,Trt,hyperparam)
  oldltarget <- logtarget(zt,Y,P,Trt,hyperparam)
  logratio <- logratio_proposal + newltarget$logtarget - oldltarget$logtarget
  
  
  
  #print("New candicate of Z: ")
  #print(as.vector(znew),quote = F)
  #decide whether accept the new value of z 
 
  #logratio <- logratio_proposal + 
  #  logtarget(znew,a0,h,c,lambda,v,w1,w2,Y,X) - 
  #  logtarget(zt,a0,h,c,lambda,v,w1,w2,Y,X)
  
  #if (!is.na(logratio) & (logratio > 0 || exp(logratio) > u)){
  if (logratio > 0 || exp(logratio) > u){
    if (verbose){cat(sprintf('Accept the new value \n'))}
    return(list(accept=TRUE,z=znew,logtargetval = newltarget$logtarget, logratio = logratio,u=u,post_beta = newltarget$post_beta))
  } 
  if (verbose) {cat(sprintf('Reject the new value \n'))}
  return(list(accept=FALSE,z=zt, logtargetval = oldltarget$logtarget,logratio = logratio,u=u,post_beta = oldltarget$post_beta))
}
