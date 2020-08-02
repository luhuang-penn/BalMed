#This is the function to perform one update iteration for all parameters
#Input are all the unkown parameters, data(X,Y,T), hyperparameters are wrapped in a list
#Return the updated value for all unkown parameters
#library(stats)
#library(MASS)
#source('../balance.R')
#source('../iteration_P.R')
#source('../iteration_z.R')


iteration <- function(X,Y,Trt,mu0,mu1,invSigma,P,z,hyperparam){
  #X is the count matrix
  alrP <- alr(P)
  invSigma_update <- itertion_invSigma(alrP,Trt,mu0,mu1,hyperparam)
  mu_update <- iteration_mu(alrP,Trt,invSigma,hyperparam)
  mu0_update <- matrix(mu_update$mu0_update)
  mu1_update <-matrix(mu_update$mu1_update)
  z_update <-iteration_z(Y,P,Trt,z,hyperparam)
  
  P_update <- iteration_P(X,Y,Trt,P,alrP,invSigma_update,mu0_update,mu1_update,z_update$z,hyperparam)
  return(list(invSigma=invSigma_update,
              mu = mu_update,
              #mu0=mu0_update,
              #mu1=mu1_update,
              z=z_update,
              P=P_update))
}
itertion_invSigma <- function(alrP,Trt,mu0,mu1,hyperparam){
  n0 <- sum(Trt==0)
  n1 <- sum(Trt==1)
  mat <- Reduce('+',lapply(1:(n0+n1),function(i) {
          mu <- ifelse(Trt[i],mu1,mu0)
          vec <- alrP[i,,drop=F]-mu
          t(vec)%*%vec}))
  mat <- solve(hyperparam$Psi+mat)
  invSigma_update <- rWishart(1, n0+n1+hyperparam$rho, mat)
  return(matrix(invSigma_update,nrow=dim(alrP)[2]))
}

iteration_mu <- function(alrP,Trt,invSigma,hyperparam){
  n0 <- sum(Trt==0)
  n1 <- sum(Trt==1)
  sum <- apply(alrP,2,function(x){tapply(x,Trt,sum)})
  sumofn1 <- sum[2,]
  varn1 <-solve(hyperparam$invOmega+n1*invSigma)
  meann1 <-varn1%*%(invSigma%*%sumofn1+hyperparam$invOmega%*%hyperparam$eta)
  mu1_update <- mvrnorm(n = 1, meann1, varn1)
  sumofn0 <- sum[1,]
  varn0 <-solve(hyperparam$invOmega+n0*invSigma)
  meann0 <-varn0%*%(invSigma%*%sumofn0+hyperparam$invOmega%*%hyperparam$eta)
  mu0_update <- mvrnorm(1,meann0,varn0)
  return(list(mu1_update=mu1_update,
              mu0_update=mu0_update,
              mu1_mean=meann1,
              mu0_mean=meann0,
              mu1_var=varn1,
              mu0_var=varn0))
}

#iteration_z<- function(Y,P,z,hyperparam){
  
#}