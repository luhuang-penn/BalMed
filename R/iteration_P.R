#returns the logtarget value from the balance equation
logtarget_zipi <- function(z,Y,P,Trt,hyperparam){
  
  B <- balance(P,z)
  
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
    (v + n) / 2 * log(b_star/2) 
  return(logtarget)
}
logtarget_Pi <- function(xi,yi,ti,pi,alrpi,mu0,mu1,invSigma,z){
 
  if (ti){mui <- mu1} else {mui <- mu0}
  logp <- sum(xi*log(pi))-sum(log(pi))-
             0.5*t(t(alrpi)-mui)%*%invSigma%*%(t(alrpi)-mui)+
              logtarget_zipi(z,yi,pi,ti,hyperparam)
      
  return(logp)
}
iteration_P <- function(X,Y,Trt,P,alrP,invSigma,
                        mu0,mu1,z,hyperparam){
  #print(list(invSigma,mu0,mu1,z,hyperparam))
  n <- dim(X)[1]
  #Y <- Y - mean(Y)
  new_P <- matrix(NA,nrow=dim(P)[1],ncol=dim(P)[2])
  accept <- vector('logical',dim(P)[1])
  logratio <- vector('numeric',dim(P)[1])
  u <-  vector('numeric',dim(P)[1])
  for (i in 1:n){
    alrpi <- alrP[i,,drop=FALSE]
    pi <- P[i,,drop=FALSE]
    xi <- X[i,,drop=FALSE]
    yi <- Y[i]
    ti <- Trt[i]
    oldlogtarget <- logtarget_Pi(xi,yi,ti,pi,alrpi,mu0,mu1,invSigma,z)
    newalrpi <- matrix(MASS::mvrnorm(n=1,mu=alrpi,Sigma = hyperparam$Sig_alrP_proprosal),nrow=1)
    newlogtarget <- logtarget_Pi(xi,yi,ti,alr_inv(newalrpi),newalrpi,mu0,mu1,invSigma,z)
    logB <-min(0, -1*mvtnorm::dmvnorm(newalrpi,mean=alrpi,sigma = hyperparam$Sig_alrP_proprosal,log=TRUE)-
           oldlogtarget+newlogtarget +
            mvtnorm::dmvnorm(alrpi,mean=newalrpi,sigma = hyperparam$Sig_alrP_proprosal,log=TRUE))
    #while (is.na(logB)){
    #  newalrpi <- matrix(MASS::mvrnorm(n=1,mu=alrpi,Sigma = hyperparam$Sig_alrP_proprosal),nrow=1)
    #  newlogtarget <- logtarget_Pi(xi,yi,ti,alr_inv(newalrpi),newalrpi,mu0,mu1,invSigma,z)
      
    #  logB <-min(0, -1*mvtnorm::dmvnorm(newalrpi,mean=alrpi,sigma = hyperparam$Sig_alrP_proprosal,log=TRUE)-
    #               oldlogtarget+newlogtarget +
    #               mvtnorm::dmvnorm(alrpi,mean=newalrpi,sigma = hyperparam$Sig_alrP_proprosal,log=TRUE))
      
    #}
    logratio[i] <- logB
    u[i] <- runif(1)
    if (!is.na(logB) && exp(logB) >= u[i]){
      new_P[i,] <- alr_inv(newalrpi)
      accept[i] <- TRUE
    }
    else{new_P[i,] <- t(pi)
    accept[i] <- FALSE}
    
  }
  return(list(new_P=new_P,accept=accept,logratio=logratio,u=u))
}
