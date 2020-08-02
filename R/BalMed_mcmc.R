#this function runs the sampling
#source("iteration.R")
#' This is the main function of balance mediation analysis
#' Runs a MCMC chain to sampling all the unknown parameters in the model
#' @param X compositional matrix
#' @param Y outcome vector
#' @param Trt treatment vector currently only accept binary value
#' @param hyperparam a list containing all the prior distribution informatiion
#' @param mu0_inital starting value for the mean in untreated group
#' @param mu1_inital starting value for the mean in the treated group
#' @param invSigma_inital starting value for the variance-covariance matrix
#' @param P_inital starting value for the latent relative abundance matrix
#' @param z_inital starting value for the balance configuration vector z
#' @param niters total number of iterations
#' @param persaveiters save chain information for each of the persaveiters iterations
#' @param savepath the path to save the current chain info
#' @param filename the name of the file to save the current chain info
#' @param logfile name of the file to store some log information
#' 
#' @expose
BalMed_mcmc <- function(X,Y,Trt,hyperparam,mu0_initial,mu1_initial,
                         invSigma_initial,P_initial,z_initial,
                         niters,persaveiters,savepath,filename,logfile='MCMClog.txt'){
  n <- dim(X)[1]
  q <- dim(X)[2]
  
  mu0 <- mu0_initial
  mu1 <- mu1_initial
  invSigma <- invSigma_initial
  P <- P_initial
  z <- z_initial
  
  mu0_chain <- matrix(NA,nrow=persaveiters,ncol=q-1)
  mu1_chain <- matrix(NA,nrow=persaveiters,ncol=q-1)
  Sigma_chain <- array(NA,c(q-1,q-1,persaveiters))
  z_chain <- matrix(NA,nrow=persaveiters,ncol=q)
  z_accept <- vector('logical',persaveiters)
  P_chain <- array(NA,c(n,q,persaveiters))
  P_accept <- matrix(NA,nrow=persaveiters,ncol=n)
  post_beta_mean <- matrix(NA,nrow = persaveiters,ncol = 3)
  post_beta_var <- array(NA,c(3,3,persaveiters))
  est_beta <- matrix(NA,nrow=persaveiters,ncol=3)
  est_bal <- matrix(NA,nrow=persaveiters,ncol=1)
  cnt <- 1
  file_cnt <- 1
  for (i in 1:niters){
    if ( (i %%1e3 == 2)) {write.table(paste0("Performing Iteration # ", i, " at ", Sys.time()),logfile,append=T)}
    
    
    update <- iteration(X,Y,Trt,mu0,mu1,invSigma,P,z,hyperparam)
    #update the unkowns for next iteration
    mu0 <- update$mu$mu0_update
    mu1 <- update$mu$mu1_update
    invSigma <- update$invSigma
    
    z <- update$z$z
    #print(paste0('updating z ',i,'th iter'))
    P <- update$P$new_P
    # record chain results for parameters
    mu0_chain[cnt,] <- mu0
    mu1_chain[cnt,] <- mu1
    Sigma_chain[,,cnt] <- solve(invSigma)
    z_chain[cnt,] <- z
    z_accept[cnt] <-update$z$accept
    P_chain[,,cnt] <-P
    P_accept[cnt,] <-update$P$accept
    # record conditional posterior dist for integrated out coefficients
    post_beta_mean[cnt,] <- update$z$post_beta$mean
    post_beta_var[,,cnt] <- update$z$post_beta$sigma
    # record estimated beta value
    bal <- balance(P,z)
    est_beta[cnt,] <- as.vector(lm(Y~Trt+bal)$coefficient)
    # record estimated a value
    est_bal[cnt,] <- balance(alr_inv(matrix(mu1-mu0,nrow=1)),z)
    cnt = cnt + 1
    if (cnt == persaveiters + 1){
      tosave <- list( mu0_chain = mu0_chain,mu1_chain = mu1_chain,
                      Sigma_chain = Sigma_chain, z_chain = z_chain,
                      z_accept = z_accept,P_chain = P_chain,P_accept = P_accept,
                      post_beta_mean = post_beta_mean,post_beta_var = post_beta_var,
                      est_beta = est_beta,est_bal=est_bal)
      save(tosave,file = paste0(savepath,filename,'_res_',file_cnt,'.RData'))
      #reinitialize
      mu0_chain <- matrix(NA,nrow=persaveiters,ncol=q-1)
      mu1_chain <- matrix(NA,nrow=persaveiters,ncol=q-1)
      Sigma_chain <- array(NA,c(q-1,q-1,persaveiters))
      z_chain <- matrix(NA,nrow=persaveiters,ncol=q)
      z_accept <- vector('logical',persaveiters)
      P_chain <- array(NA,c(n,q,persaveiters))
      P_accept <- matrix(NA,nrow=persaveiters,ncol=n)
      post_beta_mean <- matrix(NA,nrow = persaveiters,ncol = 3)
      post_beta_var <- array(NA,c(3,3,persaveiters))
      est_beta <- matrix(NA,nrow = persaveiters,ncol = 3)
      est_bal <- matrix(NA,nrow = persaveiters,ncol = 1)
      cnt <- 1
      #file_cnt increment
      file_cnt <- file_cnt + 1
    }
  }
  if (niters %% persaveiters > 0){
    tosave <- list( mu0_chain = mu0_chain,mu1_chain = mu1_chain,
                    Sigma_chain = Sigma_chain, z_chain = z_chain,
                    z_accept = z_accept,P_chain = P_chain,P_accept = P_accept,
                    post_beta_mean = post_beta_mean,post_beta_var = post_beta_var,
                    est_beta = est_beta,est_bal=est_bal)
  save(tosave,file = paste0(savepath,filename,'_res_',file_cnt,'.RData'))
  }
}

