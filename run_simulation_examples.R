#setwd("Z:/plasso/code_github/")
#.libPaths(c(.libPaths(),"Z:/R_packages"))

source("dgp.R")
source("Method.R")
source("simulation_funcs.R")

library(foreach)
library(doParallel)
library(abind)

#### basic setting ####
other_setting = read.csv("other_setting3.csv")
repeat_times = 100

#### ####
cl=makeCluster(min(20, repeat_times), outfile="")
registerDoParallel(cl)
set.seed(20201110)
seed_each = sample.int(50000, repeat_times)

#### ####
for(example in 1:5){
  for(j in 1:4){                 # determin n,p and rho
    start.time = Sys.time()
    print(start.time)
    
    #### initial settings ####
    n = other_setting$n[j]
    p = other_setting$p[j]
    rou = other_setting$rou[j]
    n_val=1000
    nfolds=4
    onestd = F
    
    #### settings about true beta ####
    set.seed(100)
    beta = c(-1.5,runif(15,-0.7,0.7),rep(0,p-15))
    
    print(sprintf("==== Example %d, j=%d: (n,p)=(%d,%d), rou=%.2f ====", example, j, n, p, rou))
    
    #### simulation (parallel) ####
    if(example==1|example==5){
      #### settings about beta_hat (prior information) ####
      n_prior = 3
      noisy_levels = c(0.25,0.25,0.5)
      K = 3
      beta_hat = matrix(NA,p+1,K)
      for(k in 1:K){
        noisy = noisy_levels[k]
        if(example == 1){
          beta_hat[,k] = c(beta[1],beta[2:31]+runif(30,-noisy,noisy),beta[32:(p+1)])
        }else if(example == 5){
          beta_hat[,k] = c(-1.5,runif(p,-1.5,1.5))
        }
      }
      #### computation ####
      result_all = foreach(seed=seed_each,.combine=sim.combine) %dopar% 
        simulation1(seed=seed, n=n, p=p, n_val=n_val, rou=rou,
                    beta=beta, beta_hat=beta_hat,
                    metric="loglike_val", onestd = F,nfolds = nfolds
        )
      
    }else if(example == 2){
      #### settings about beta support (prior information) ####
      n_prior = 3
      s1 = c(rep(T,7),rep(F,p-7))
      s2 = c(rep(F,4),rep(T,7),rep(F,p-11))
      s3 = c(rep(F,8),rep(T,7),rep(F,p-15))
      s1[16:17] = T
      s2[18:19] = T
      s3[20:21] = T
      beta_supports = cbind(s1,s2,s3)
      
      #### calculation ####
      result_all = foreach(seed=seed_each,.combine=sim.combine) %dopar% 
        simulation2(seed=seed, n=n, p=p, n_val=n_val, rou=rou,
                    beta=beta, beta_supports=beta_supports,
                    metric = "loglike_val", onestd = F,nfolds = nfolds
        )
      
    }else if(example == 3){
      n_prior=4
      #### settings about beta_hat and beta_support (prior information) ####
      noisy_levels = c(0.25,0.25)
      K1 = 2
      beta_hat = matrix(NA,p+1,K1)
      for(k in 1:K1){
        noisy = noisy_levels[k]
        beta_hat[,k] = c(beta[1],beta[2:31]+runif(30,-noisy,noisy),beta[32:(p+1)])
      }
      
      #### setting about beta_supports ####
      s1 = c(rep(T,5),rep(F,5),rep(T,5),rep(F,p-15))
      s2 = c(rep(F,5),rep(T,5),rep(T,5),rep(F,p-15))
      s1[16:17] = T
      s2[18:19] = T
      beta_supports = cbind(s1,s2)
      
      #### calculation ####
      result_all = foreach(seed=seed_each,.combine=sim.combine) %dopar% 
        simulation3(seed=seed,n=n,p=p,n_val=n_val,rou=rou,
                    beta=beta,beta_hat=beta_hat,beta_supports=beta_supports,
                    metric="loglike_val",onestd = F,nfolds = nfolds
                    )
      
    }else if(example == 4){
      n_prior = 3
      #### prior information are extracted from prior datasets ####
      #### calculation ####
      result_all = foreach(seed=seed_each,.combine=sim.combine) %dopar% 
        simulation4(seed=seed, n=n, p=p, n_val=n_val, rou=rou,
                    beta=beta, n_prior_samples=c(2000,3000,4000), h=0.2,
                    metric="loglike_val",onestd = F, nfolds = nfolds)
    }
    
    
    #### output ####
    result_mean = apply(result_all$result_out,2:3,mean,na.rm=T)
    result_std = apply(result_all$result_out,2:3,sd,na.rm=T)
    
    result_out0 = result.table(result_mean,result_std)
    result_out1 = data.frame(result_out0,
                             row.names=c("Original",paste("prior",1:n_prior,sep=""),"Our Method")
    )
    colnames(result_out1) <- c("FNR","FPR","FR","CNZ","INZ","MCC","FDR",
                               "BETA_MSE","BETA_BIAS",
                               "ACCURACY_val","AUC_val","loglike_val")
    
    write.csv(result_out1,file=sprintf("./result/example%d_%d_result.csv", example, j))
    write.csv(result_all$tuning,file=sprintf("./result/example%d_%d_tunings.csv", example, j))
    
    ####
    end.time = Sys.time()
    print(end.time)
    print(end.time-start.time)
    
  }
}



stopCluster(cl)





