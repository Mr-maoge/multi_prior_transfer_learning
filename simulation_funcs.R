source("dgp.R")
source("Method.R")
library(abind)


#### part 1: prior information is beta coefficient ####
simulation1 <- function(seed=NULL,n=1000, p=100, n_val=1000, rou=0.5,
                        beta=NULL, beta_hat=NULL,
                        metric="loglike_val",onestd = F,nfolds = 5)
{
  #### Function for example 1 and 5 in paper, prior information is the beta coefficients ####

  # Args:
  #  seed: random seed
  #  n: sample size 
  #  p: num of variables
  #  n_val: sample size for validation set 
  #  rou: parameter for the AR structure predictor covariance matrix 
  #  beta: true coefficient (p, )
  #  beta_hat: prior beta information (p, K)  
  #  metric: metric for CV to select the best model, one of "loglike_val","AUC_val"
  #  onestd: logic, whether "mean + 1 std" strategy is used to select the best model (default False) 
  #  nfolds: number of folds in CV  
  
  # Return:
  #  result_out: result table (1, K+2, 12)
  #              dim 2 for different methods: (Origin, prior-1, ..., prior-K, Propsed) 
  #              dim 3 for evaluation metrics: ("FNR","FPR","FR","CNZ","INZ","MCC","FDR", "BETA_MSE","BETA_BIAS", "ACCURACY_val","AUC_val","loglike_val")
  
  #  tuning: tuning parameters for different methods
  
  source("Method.R")
  source("dgp.R")
  
  lambda_best = NA
  v_best = NA
  eta = NA
  tryCatch({
    #### initialize ####
    if(is.null(beta)){
      beta=c(-1,
             -0.1,-1,2,-0.2,0.1,
             0.1,1,-2,0.2,-0.1,
             -0.1,-1,2,-0.2,0.1,rep(0,p-15))
    }
    
    if(is.null(beta_hat)){
      beta_hat_1 = beta + runif(p+1,-0.01,0.01)
      beta_hat_2 = beta + runif(p+1,-0.01,0.01)
      beta_hat_3 = beta + runif(p+1,-0.01,0.01)
      beta_hat = cbind(beta_hat_1,beta_hat_2,beta_hat_3)
    }
    
    penalty = "SCAD"
    # lambda_seq : candidate lambdas for methods except our method
    lambda_seq = round(exp(seq(log(0.2), log(0.01),length.out=10)),digits=5)
    # lambda_seq : candidate lambdas for our method
    lambda_seq2 = round(exp(seq(log(0.2), log(0.01),length.out=10)),digits=5)
    # lambda_seq : candidate v for our method
    a=seq(0,0.9,0.15)
    v_seq = v = a/(1-a)
    
    if(!is.null(seed)){
      set.seed(seed)
    }
    
    n_prior = dim(beta_hat)[2]
    result = matrix(NA,5,12)
    result_out = array(NA,dim=c(1,5,12))
    
    #### generate data ####
    K = dim(beta_hat)[2]
    data = data_generator1(n,p,n_val=n_val,rou = rou,beta=beta)
    x = data$x
    y = data$y
    x_val = data$x_val
    y_val = data$y_val
    beta.oracle = data$beta
    print("Generate data end")
    
    #### create fold ####
    idx0 = which(y<0.5)
    idx1 = which(y>0.5)
    tmp0 = createFolds(1:length(idx0),k=nfolds,list=TRUE,returnTrain=F)
    tmp1 = createFolds(1:length(idx1),k=nfolds,list=TRUE,returnTrain=F)
    fold = rep(NA,n)
    for(i in 1:nfolds){
      idx = c(idx0[tmp0[[i]]],idx1[tmp1[[i]]])
      fold[idx] = i
    }
    
    #### train models ####
    #### Don't use prior infomation ####
    m0 = cv.logistic(x,y,penalty = penalty,
                     lambda=lambda_seq,nfolds = nfolds,fold=fold,
                     metric=metric,onestd=onestd,
                     x_val=x_val,y_val=y_val,beta.oracle = beta.oracle)
    m0_best = m0$best
    m0_best_lambda = m0_best$lambda
    m0_best_metric = m0_best$metric_val
    cat("lambda_best=",m0_best_lambda,"metric=",m0_best_metric,"\n")
    print("m0 end")
    
    #### sigle source prior transfer method ####
    y_hat = matrix(NA,n,K)
    m1_loglike = rep(NA,K)
    m1_best = list()
    for(k in 1:K){
      beta_hat_tmp = beta_hat[,k]
      y_hat[,k] = expit(beta_hat_tmp[1]+x%*%beta_hat_tmp[-1])[,1]
      m1_loglike[k] = loglike.logistic(y,y_hat[,k])
      
      m_tmp2 = cv.transfer.logistic(x, y, matrix(y_hat[,k],length(y),1), v_seq, eta0=1,
                                    penalty=penalty, lambda=lambda_seq2, nfolds = nfolds,fold=fold,
                                    metric=metric,onestd = onestd,
                                    x_val=x_val,y_val=y_val,beta.oracle = beta.oracle)
      m_tmp2_best = m_tmp2$best
      m1_best[[k]] = m_tmp2_best
      cat("lambda_best=",m_tmp2_best$lambda,"metric=",m_tmp2_best$metric_val,"\n")
    }
    print("m1 end")
    
    #### our model ####
    eta0 = weight.prior(m1_loglike,1)
    m2 = cv.transfer.logistic(x,y,y_hat,v_seq,eta0,penalty,lambda_seq2,
                              nfolds = nfolds, fold=fold,
                              metric=metric,onestd = onestd,
                              x_val=x_val,y_val=y_val,beta.oracle = beta.oracle)
    m2_best = m2$best
    eta = m2_best$eta
    cat("lambda_best=",m2_best$lambda,"metric=",m2_best$metric_val,"\n")
    print(eta)
    print("m2 end")
    
    #### output ####
    m_best_sets = list()
    m_best_sets[[1]] = m0_best
    for(i in 1:n_prior){
      m_best_sets[[i+1]] = m1_best[[i]]
    }
    m_best_sets[[n_prior+2]] = m2_best
    tuning = c()
    for(i in 1:(2+n_prior)){
      m_best_tmp = m_best_sets[[i]]
      result[i,] = c(m_best_tmp$FNR,m_best_tmp$FPR,m_best_tmp$FR,
                     m_best_tmp$CNZ,m_best_tmp$INZ,m_best_tmp$MCC,m_best_tmp$FDR,
                     m_best_tmp$BETA_MSE,m_best_tmp$BETA_BIAS,
                     m_best_tmp$ACURACCY_val,m_best_tmp$AUC_val,m_best_tmp$loglike_val)
      tuning = c(tuning,m_best_tmp$lambda,m_best_tmp$eta)
    }
    #### end ####
    result_out[1,,] = result
    return(list(result_out=result_out,tuning=tuning))
  },error=function(e){
    print(e$message)
  })
  return(list(result_out=NULL,tuning=NULL)) # if error occurs
}

#### part 2: prior information is beta_support ####
simulation2 <- function(seed=NULL,n=1000, p=100, n_val=1000, rou=0.5, 
                        beta=NULL, beta_supports=NULL,
                        metric = "loglike_val",onestd = F,nfolds = 5)
{
  #### Function for example 2 in paper, prior information is the beta supports ####
  
  # Args:
  #  seed: random seed
  #  n: sample size 
  #  p: num of variables
  #  n_val: sample size for validation set 
  #  rou: parameter for the AR structure predictor covariance matrix 
  #  beta: true coefficient (p, )
  #  beta_supports: prior beta information (p, K), `1` or `T` in the matrix indicates the prior information show that the related variable is important 
  #  metric: metric for CV to select the best model, one of "loglike_val","AUC_val"
  #  onestd: logic, whether "mean + 1 std" strategy is used to select the best model (default False) 
  #  nfolds: number of folds in CV  
  
  # Return:
  #  result_out: result table (1, K+2, 12)
  #              dim 2 for different methods: (Origin, prior-1, ..., prior-K, Propsed) 
  #              dim 3 for evaluation metrics: ("FNR","FPR","FR","CNZ","INZ","MCC","FDR", "BETA_MSE","BETA_BIAS", "ACCURACY_val","AUC_val","loglike_val")
  
  #  tuning: tuning parameters for different methods
  
  source("dgp.R")
  source("Method.R")
  
  lambda_best = NA
  v_best = NA
  eta = NA
  tryCatch({
  #### initialize ####
  if(is.null(beta)){
    beta=c(0,
           -0.5,-0.2,0.3,-0.2,0.4,
           0.5,1,-0.3,0.2,-0.4,
           -0.4,-0.5,0.3,-0.5,0.4,rep(0,p-15))
  }
    
  if(is.null(beta_supports)){
    s1 = c(rep(T,5),rep(F,5),rep(F,5),rep(F,2),rep(F,2),rep(F,p-19))
    s2 = c(rep(F,5),rep(T,5),rep(F,5),rep(F,2),rep(F,2),rep(F,p-19))
    s3 = c(rep(F,5),rep(F,5),rep(T,5),rep(F,2),rep(F,2),rep(F,p-19))
    beta_supports = cbind(s1,s2,s3)
  }  
    
  penalty = "SCAD"
  # lambda_seq : candidate lambdas for methods except our method
  lambda_seq = round(exp(seq(log(0.2), log(0.02),length.out=6)),digits=5)
  # lambda_seq : candidate lambdas for our method
  lambda_seq2 = round(exp(seq(log(0.2), log(0.02),length.out=6)),digits=5)
  # lambda_seq : candidate v for our method
  a=seq(0,0.9,0.15)
  v_seq = v = a/(1-a)
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  n_prior = dim(beta_supports)[2]
  result = matrix(NA,2+n_prior,12)
  result_out = array(NA,dim=c(1,2+n_prior,12))
  
  #### generate data ####
  K = dim(beta_supports)[2]
  data = data_generator1(n,p,n_val=n_val,n_test = 1, rou=rou, beta=beta)
  x = data$x
  y = data$y
  x_val = data$x_val
  y_val = data$y_val
  beta.oracle = data$beta
  print("Generate data end")
  
  #### create fold ####
  idx0 = which(y<0.5)
  idx1 = which(y>0.5)
  tmp0 = createFolds(1:length(idx0),k=nfolds,list=TRUE,returnTrain=F)
  tmp1 = createFolds(1:length(idx1),k=nfolds,list=TRUE,returnTrain=F)
  fold = rep(NA,n)
  for(i in 1:nfolds){
    idx = c(idx0[tmp0[[i]]],idx1[tmp1[[i]]])
    fold[idx] = i
  }
  
  #### train models ####
  #### Don't use prior infomation ####
  m0 = cv.logistic(x,y,penalty = penalty,
                   lambda=lambda_seq,nfolds = nfolds,fold=fold,
                   metric=metric,onestd=onestd,
                   x_val=x_val,y_val=y_val,beta.oracle = beta.oracle)
  m0_best = m0$best
  m0_best_lambda = m0_best$lambda
  m0_best_metric = m0_best$metric_val
  cat("lambda_best=",m0_best_lambda,"metric=",m0_best_metric,"\n")
  print("m0 end")
  
  #### Fully trust prior infomation & prior-k (also construct y_hat) ####
  y_hat = matrix(NA,n,K)
  m1_best_lambda = rep(NA,K)
  m1_cv_loglike_val = rep(NA,K)
  m1_best = list()
  for(k in 1:K){
    beta_support = beta_supports[,k]
    m_tmp = cv.logistic(x,y,penalty = penalty,beta_support = beta_support,
                        lambda=lambda_seq,nfolds = nfolds,fold=fold,
                        metric=metric,onestd=onestd,
                        x_val=x_val,y_val=y_val,beta.oracle = beta.oracle)
    m1_cv_loglike_val[k] = m_tmp$loglike_val_cv
    m_tmp_best = m_tmp$best
    beta = as.vector(m_tmp_best$beta)
    y_hat_tmp = as.vector(expit(x%*%beta[-1]+beta[1]))
    y_hat[,k] = y_hat_tmp
    
    # prior lasso method (1 source)
    m_tmp2 = cv.transfer.logistic(x,y,matrix(y_hat[,k],length(y),1),v_seq,eta=1,
                                  penalty,lambda_seq2,nfolds = nfolds,fold=fold,
                                  metric=metric,onestd = onestd,
                                  x_val=x_val,y_val=y_val,beta.oracle = beta.oracle)
    m_tmp2_best = m_tmp2$best
    m1_best[[k]] = m_tmp2_best
    cat("lambda_best=",m_tmp2_best$lambda,"metric=",m_tmp2_best$metric_val,"\n")
  }
  print("m1 end")
  
  #### our model ####
  eta0 = weight.prior(m1_cv_loglike_val,1)
  m2 = cv.transfer.logistic(x,y,y_hat,v_seq,eta0,penalty,lambda_seq2,nfolds = nfolds,fold=fold,
                            metric=metric,onestd = onestd,
                            x_val=x_val,y_val=y_val,beta.oracle = beta.oracle)
  m2_best = m2$best
  eta = m2_best$eta
  cat("lambda_best=",m2_best$lambda,"metric=",m2_best$metric_val,"\n")
  print(eta)
  print("m2 end")
  
  #### output ####
  m_best_sets = list()
  m_best_sets[[1]] = m0_best
  for(i in 1:n_prior){
    m_best_sets[[i+1]] = m1_best[[i]]
  }
  m_best_sets[[n_prior+2]] = m2_best
  tuning = c()
  for(i in 1:(2+n_prior)){
    m_best_tmp = m_best_sets[[i]]
    result[i,] = c(m_best_tmp$FNR,m_best_tmp$FPR,m_best_tmp$FR,
                   m_best_tmp$CNZ,m_best_tmp$INZ,m_best_tmp$MCC,m_best_tmp$FDR,
                   m_best_tmp$BETA_MSE,m_best_tmp$BETA_BIAS,
                   m_best_tmp$ACURACCY_val,m_best_tmp$AUC_val,m_best_tmp$loglike_val)
    tuning = c(tuning,m_best_tmp$lambda,m_best_tmp$eta)
  }
  #### end ####
  result_out[1,,] = result
  return(list(result_out=result_out,
              tuning=tuning))
  },error=function(e){
    print(e$message)
  })
  return(list(result_out=NULL,tuning=NULL))
}

#### part 3: 2 beta_hat + 2 beta support
simulation3 <- function(seed=NULL, n=1000, p=100, n_val=1000, rou=0.5,
                        beta=NULL, beta_hat=NULL, beta_supports=NULL,
                        metric="loglike_val",onestd = F,nfolds = 5)
{
  #### Function for example 3 in paper, prior information is the beta hat + beta supports ####
  
  # Args:
  #  seed: random seed
  #  n: sample size 
  #  p: num of variables
  #  n_val: sample size for validation set 
  #  rou: parameter for the AR structure predictor covariance matrix 
  #  beta: true coefficient (p, )
  #  beta_hat: prior beta information (p, K1)  
  #  beta_supports: prior beta information (p, K2), `1` or `T` in the matrix indicates the prior information show that the related variable is important 
  #  metric: metric for CV to select the best model, one of "loglike_val","AUC_val"
  #  onestd: logic, whether "mean + 1 std" strategy is used to select the best model (default False) 
  #  nfolds: number of folds in CV  
  
  # Return:
  #  result_out: result table (1, K+2, 12)
  #              dim 2 for different methods: (Origin, prior-1, ..., prior-K, Propsed) 
  #              dim 3 for evaluation metrics: ("FNR","FPR","FR","CNZ","INZ","MCC","FDR", "BETA_MSE","BETA_BIAS", "ACCURACY_val","AUC_val","loglike_val")
  
  #  tuning: tuning parameters for different methods
  
  source("Method.R")
  source("dgp.R")
  
  lambda_best = NA
  v_best = NA
  eta = NA
  tryCatch({
    #### initialize ####
    if(is.null(beta)){
      beta=c(-1,
             -0.1,-1,2,-0.2,0.1,
             0.1,1,-2,0.2,-0.1,
             -0.1,-1,2,-0.2,0.1,rep(0,p-15))
    }
    
    if(is.null(beta_hat)){
      beta_hat_1 = beta + runif(p+1,-0.01,0.01)
      beta_hat_2 = beta + runif(p+1,-0.01,0.01)
      beta_hat = cbind(beta_hat_1,beta_hat_2)
    }
    
    if(is.null(beta_supports)){
      s1 = c(rep(T,5),rep(F,5),rep(F,5),rep(F,2),rep(F,2),rep(F,p-19))
      s2 = c(rep(F,5),rep(T,5),rep(F,5),rep(F,2),rep(F,2),rep(F,p-19))
      beta_supports = cbind(s1,s2)
    }  
    
    penalty = "SCAD"
    # lambda_seq : candidate lambdas for methods except our method
    lambda_seq = round(exp(seq(log(0.2), log(0.01),length.out=10)),digits=5)
    # lambda_seq : candidate lambdas for our method
    lambda_seq2 = round(exp(seq(log(0.2), log(0.01),length.out=10)),digits=5)
    # lambda_seq : candidate v for our method
    a=seq(0,0.9,0.15)
    v_seq = v = a/(1-a)
    
    if(!is.null(seed)){
      set.seed(seed)
    }
    
    K1 = dim(beta_hat)[2]
    K2 = dim(beta_supports)[2]
    K = K1+K2
    
    result = matrix(NA,2+K,12)
    result_out = array(NA,dim=c(1,2+K,12))
    
    #### generate data ####
    data = data_generator1(n,p,n_val=n_val,rou = rou,beta=beta)
    x = data$x
    y = data$y
    x_val = data$x_val
    y_val = data$y_val
    beta.oracle = data$beta
    print("Generate data end")
    
    #### create fold ####
    idx0 = which(y<0.5)
    idx1 = which(y>0.5)
    tmp0 = createFolds(1:length(idx0),k=nfolds,list=TRUE,returnTrain=F)
    tmp1 = createFolds(1:length(idx1),k=nfolds,list=TRUE,returnTrain=F)
    fold = rep(NA,n)
    for(i in 1:nfolds){
      idx = c(idx0[tmp0[[i]]],idx1[tmp1[[i]]])
      fold[idx] = i
    }
    
    #### train models ####
    #### Don't use prior infomation ####
    m0 = cv.logistic(x,y,penalty = penalty,
                     lambda=lambda_seq,nfolds = nfolds,fold=fold,
                     metric=metric,onestd=onestd,
                     x_val=x_val,y_val=y_val,beta.oracle = beta.oracle)
    m0_best = m0$best
    m0_best_lambda = m0_best$lambda
    m0_best_metric = m0_best$metric_val
    cat("lambda_best=",m0_best_lambda,"metric=",m0_best_metric,"\n")
    print("m0 end")
    #### Fully trust prior infomation (also construct y_hat) ####
    y_hat = matrix(NA,n,K)
    m1_loglike = rep(NA,K)
    m1_best = list()
    # beta hat part
    for(k in 1:K1){
      beta_hat_tmp = beta_hat[,k]
      y_hat[,k] = expit(beta_hat_tmp[1]+x%*%beta_hat_tmp[-1])[,1]
      m1_loglike[k] = loglike.logistic(y,y_hat[,k])
      
      m_tmp2 = cv.transfer.logistic(x, y, matrix(y_hat[,k],length(y),1), v_seq, eta=1,
                                    penalty,lambda_seq2,nfolds = nfolds,fold=fold,
                                    metric=metric,onestd = onestd,
                                    x_val=x_val,y_val=y_val,beta.oracle = beta.oracle)
      m_tmp2_best = m_tmp2$best
      m1_best[[k]] = m_tmp2_best
      cat("lambda_best=",m_tmp2_best$lambda,"metric=",m_tmp2_best$metric_val,"\n")
    }
    print("m1 1 end") 
    
    # beta supports part
    for(k in 1:K2){
      beta_support = beta_supports[,k]
      m_tmp = cv.logistic(x,y,penalty = penalty,beta_support = beta_support,
                          lambda=lambda_seq,nfolds = nfolds,fold=fold,
                          metric=metric,onestd=onestd,
                          x_val=x_val,y_val=y_val,beta.oracle = beta.oracle)
      m_tmp_best = m_tmp$best
      beta = as.vector(m_tmp_best$beta)
      y_hat_tmp = as.vector(expit(x%*%beta[-1]+beta[1]))
      y_hat[,k+K1] = y_hat_tmp
      
      m1_loglike[k+K1] = m_tmp$loglike_val_cv
      # prior lasso method (1 source)
      m_tmp2 = cv.transfer.logistic(x,y,matrix(y_hat[,k+K1],length(y),1),v_seq,eta=1,
                                    penalty,lambda_seq2,nfolds = nfolds,fold=fold,
                                    metric=metric,onestd = onestd,
                                    x_val=x_val,y_val=y_val,beta.oracle = beta.oracle)
      m_tmp2_best = m_tmp2$best
      m1_best[[k+K1]] = m_tmp2_best
      cat("lambda_best=",m_tmp2_best$lambda,"metric=",m_tmp2_best$metric_val,"\n")
    }
    print("m1 2 end")
    
    #### our model ####
    eta0 = weight.prior(m1_loglike,1)
    m2 = cv.transfer.logistic(x,y,y_hat,v_seq, eta0, penalty, lambda_seq2, nfolds = nfolds,fold=fold,
                              metric=metric,onestd = onestd,
                              x_val=x_val,y_val=y_val,beta.oracle = beta.oracle)
    m2_best = m2$best
    eta = m2_best$eta
    cat("lambda_best=",m2_best$lambda,"metric=",m2_best$metric_val,"\n")
    print(eta)
    print("m2 end")
    
    #### output ####
    m_best_sets = list()
    m_best_sets[[1]] = m0_best
    for(i in 1:K){
      m_best_sets[[i+1]] = m1_best[[i]]
    }
    m_best_sets[[K+2]] = m2_best
    tuning = c()
    for(i in 1:(2+K)){
      m_best_tmp = m_best_sets[[i]]
      result[i,] = c(m_best_tmp$FNR,m_best_tmp$FPR,m_best_tmp$FR,
                     m_best_tmp$CNZ,m_best_tmp$INZ,m_best_tmp$MCC,m_best_tmp$FDR,
                     m_best_tmp$BETA_MSE,m_best_tmp$BETA_BIAS,
                     m_best_tmp$ACURACCY_val,m_best_tmp$AUC_val,m_best_tmp$loglike_val)
      tuning = c(tuning,m_best_tmp$lambda,m_best_tmp$eta)
    }
    #### end ####
    result_out[1,,] = result
    return(list(result_out=result_out,tuning=tuning))
  },error=function(e){
    print(e$message)
  })
  return(list(result_out=NULL,tuning=NULL)) # if error occurs
}


#### part 4: prior information is derived from dataset (beta coeffcients) ####
simulation4 <- function(seed=NULL,n=1000,p=100,n_val=1000,rou=0.5,
                        beta=NULL,n_prior_samples=c(2000,2000,2000),h=0.1,
                        metric="loglike_val",onestd = F,nfolds = 5)
{
  #### Function for example 3 in paper, prior information is the beta hat + beta supports ####
  
  # Args:
  #  seed: random seed
  #  n: sample size 
  #  p: num of variables
  #  n_val: sample size for validation set 
  #  rou: parameter for the AR structure predictor covariance matrix 
  #  beta: true coefficient (p, )
  
  #  n_prior_samples: sample size for prior datasets   
  #  h: noise level for prior 
  
  #  metric: metric for CV to select the best model, one of "loglike_val","AUC_val"
  #  onestd: logic, whether "mean + 1 std" strategy is used to select the best model (default False) 
  #  nfolds: number of folds in CV  
  
  # Return:
  #  result_out: result table (1, K+2, 12)
  #              dim 2 for different methods: (Origin, prior-1, ..., prior-K, Propsed) 
  #              dim 3 for evaluation metrics: ("FNR","FPR","FR","CNZ","INZ","MCC","FDR", "BETA_MSE","BETA_BIAS", "ACCURACY_val","AUC_val","loglike_val")
  
  #  tuning: tuning parameters for different methods
  
  source("Method.R")
  source("dgp.R")
  
  lambda_best = NA
  v_best = NA
  eta = NA
  tryCatch({
    #### initialize ####
    if(is.null(beta)){
      beta=c(0,
             -0.1,-1,2,-0.2,0.1,
             0.1,1,-2,0.2,-0.1,
             -0.1,-1,2,-0.2,0.1,rep(0,p-15))
    }
    
    penalty = "SCAD"
    # lambda_seq : candidate lambdas for methods except our method
    lambda_seq = round(exp(seq(log(0.2), log(0.01),length.out=10)),digits=5)
    # lambda_seq : candidate lambdas for our method
    lambda_seq2 = round(exp(seq(log(0.2), log(0.01),length.out=10)),digits=5)
    # lambda_seq : candidate v for our method
    a=seq(0,0.9,0.15)
    v_seq = v = a/(1-a)
    
    if(!is.null(seed)){
      set.seed(seed)
    }
    
    K = length(n_prior_samples)
    result = matrix(NA,2+K,12)
    result_out = array(NA,dim=c(1,2+K,12))
    
    #### generate data ####
    data = data_generator1(n,p,n_val=n_val,rou=rou,beta=beta)
    x = data$x
    y = data$y
    x_val = data$x_val
    y_val = data$y_val
    beta.oracle = data$beta
    print("Generate data end")
    
    #### create fold ####
    idx0 = which(y<0.5)
    idx1 = which(y>0.5)
    tmp0 = createFolds(1:length(idx0),k=nfolds,list=TRUE,returnTrain=F)
    tmp1 = createFolds(1:length(idx1),k=nfolds,list=TRUE,returnTrain=F)
    fold = rep(NA,n)
    for(i in 1:nfolds){
      idx = c(idx0[tmp0[[i]]],idx1[tmp1[[i]]])
      fold[idx] = i
    }
    
    #### train models ####
    #### Don't use prior infomation ####
    m0 = cv.logistic(x,y,penalty = penalty,
                     lambda=lambda_seq,nfolds = nfolds,fold=fold,
                     metric=metric,onestd=onestd,
                     x_val=x_val,y_val=y_val,beta.oracle = beta.oracle)
    m0_best = m0$best
    m0_best_lambda = m0_best$lambda
    m0_best_metric = m0_best$metric_val
    cat("lambda_best=",m0_best_lambda,"metric=",m0_best_metric,"\n")
    print("m0 end")
    
    #### Fully trust prior infomation (also construct y_hat) ####
    # generate prior dataset + generate beta_hat
    beta_hat = matrix(NA,p+1,K)
    for(k in 1:K){
      beta_prior = beta + c(0,runif(30,-h,h),rep(0,p-30))
      data_prior = data_generator1(n_prior_samples[k],p,n_val=1,n_test=1,
                                   rou=rou,beta=beta_prior)
      m_tmp = cv.logistic(data_prior$x,data_prior$y,penalty = penalty,
                       lambda=lambda_seq,nfolds = nfolds,fold=fold,
                       metric=metric,onestd=onestd)$best
      beta_hat[,k] = m_tmp$beta
    }
    cat("Generate beta from datasets ends.","\n")
    ## the same procedures as if we known beta_hat in advance
    y_hat = matrix(NA,n,K)
    m1_loglike = rep(NA,K)
    m1_best = list()
    for(k in 1:K){
      beta_hat_tmp = beta_hat[,k]
      y_hat[,k] = expit(beta_hat_tmp[1]+x%*%beta_hat_tmp[-1])[,1]
      m1_loglike[k] = loglike.logistic(y,y_hat[,k])
      
      m_tmp2 = cv.transfer.logistic(x,y,matrix(y_hat[,k],length(y),1),v_seq,eta=1,
                                    penalty,lambda_seq2,nfolds = nfolds,fold=fold,
                                    metric=metric,onestd = onestd,
                                    x_val=x_val,y_val=y_val,beta.oracle = beta.oracle)
      m_tmp2_best = m_tmp2$best
      m1_best[[k]] = m_tmp2_best
      cat("lambda_best=",m_tmp2_best$lambda,"metric=",m_tmp2_best$metric_val,"\n")
    }
    print("m1 end")
    #### our model ####
    eta0 = weight.prior(m1_loglike,1)
    m2 = cv.transfer.logistic(x,y,y_hat,v_seq,eta0,penalty,lambda_seq2,nfolds = nfolds,fold=fold,
                              metric=metric,onestd = onestd,
                              x_val=x_val,y_val=y_val,beta.oracle = beta.oracle)
    m2_best = m2$best
    eta = m2_best$eta
    cat("lambda_best=",m2_best$lambda,"metric=",m2_best$metric_val,"\n")
    print(eta)
    print("m2 end")
    
    #### output ####
    m_best_sets = list()
    m_best_sets[[1]] = m0_best
    for(i in 1:K){
      m_best_sets[[i+1]] = m1_best[[i]]
    }
      m_best_sets[[K+2]] = m2_best
    tuning = c()
    for(i in 1:(2+K)){
      m_best_tmp = m_best_sets[[i]]
      result[i,] = c(m_best_tmp$FNR,m_best_tmp$FPR,m_best_tmp$FR,
                     m_best_tmp$CNZ,m_best_tmp$INZ,m_best_tmp$MCC,m_best_tmp$FDR,
                     m_best_tmp$BETA_MSE,m_best_tmp$BETA_BIAS,
                     m_best_tmp$ACURACCY_val,m_best_tmp$AUC_val,m_best_tmp$loglike_val)
      tuning = c(tuning,m_best_tmp$lambda,m_best_tmp$eta)
    }
    #### end ####
    result_out[1,,] = result
    return(list(result_out=result_out,tuning=tuning))
  },error=function(e){
    print(e$message)
  })
  return(list(result_out=NULL,tuning=NULL)) # if error occurs
}



#### some useful functions to summarize the results ####
sim.combine <- function(a,b)
{
  result_out = abind(a$result_out,b$result_out,along=1)
  tuning = rbind(a$tuning,b$tuning)
  return(list(result_out=result_out,tuning=tuning))
}

sim.combine3 <- function(a,b)
{
  result_out = abind(a$result_out,b$result_out,along=1)
  tuning = rbind(a$tuning,b$tuning)
  beta_supports = abind(a$beta_supports,b$beta_supports)
  return(list(result_out=result_out,tuning=tuning,beta_supports=beta_supports))
}

result.table <- function(r.mean,r.std)
{
  n = dim(r.mean)[1]
  p = dim(r.mean)[2]
  r = matrix(" ",n,p)
  for(i in 1:n){
    for(j in 1:p){
      r[i,j] = sprintf("%.3f(%.3f)",r.mean[i,j],r.std[i,j])
    }
  }
  r
}


