library(pROC)
library(caret)

#### functions used to fit penalized logistic regression ####
expit <- function(x)
{
  1/(1+exp(-x))
}

soft.thresholding <- function(x, lambda)
{
  sign(x)*max(abs(x)-lambda,0)
}

scad.thresholding <- function(x, lambda, gamma=3.7)
{
  #print(sprintf("%.4f,%.4f,%.4f",x,lambda,gamma))
  if(abs(x)<=2*lambda){
    x_up = soft.thresholding(x,lambda)
  }else if(abs(x)>gamma*lambda){
    x_up = x
  }else{
    x_up = soft.thresholding(x,gamma*lambda/(gamma-1))/(1-1/(gamma-1))
  }
  return(x_up)
}

mcp.thresholding <- function(x, lambda, gamma=3)
{
  if(abs(x)<=gamma*lambda){
    x_up = soft.thresholding(x,lambda)/(1-1/gamma)
  }else{
    x_up = x
  }
  return(x_up)
}

#### functions used for evaluation ####
loglike.logistic <- function(y,y_hat)
{
  y_hat = sapply(y_hat,function(x) max(1e-5,min(1-1e-5,x)))
  a = sum(y*log(y_hat) + (1-y)*log(1-y_hat))/length(y)
  return(a)
}

cal.coeff_metrics <- function(beta, beta_hat)
{
  # Function used to calculate the metrics regarding variable selections, estimatios etc.
  # Args:
  #  beta: the true coefficients vector
  #  beta_hat: the estimated coefficients vector
  #
  # Return:
  #  result: a list containing the evaluation metrics
  #
  # * both beta and beta_hat doesn't include the intercept term
  
  #### variable selection related ####
  TP = sum((beta!=0)&(beta_hat!=0))
  FP = sum((beta==0)&(beta_hat!=0))
  TN = sum((beta==0)&(beta_hat==0))
  FN = sum((beta!=0)&(beta_hat==0))

  FNR = FN / (TP+FN)
  FPR = FP / (TN+FP)
  FR = FNR + FPR
  
  CNZ = TP
  INZ = FP
  
  MCC = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  
  FDR = FP / (TP+FP)
  
  #### bias or mse ####
  BETA_BIAS = sum(abs(beta-beta_hat)*(beta!=0))
  BETA_MSE = sum((beta-beta_hat)^2) / length(beta)
  
  result = list(FNR=FNR,FPR=FPR,FR=FR,
                CNZ=CNZ,INZ=INZ,MCC=MCC,FDR=FDR,BETA_BIAS=BETA_BIAS,BETA_MSE=BETA_MSE)
  return(result)
}

cal.auc <- function(y,y_hat)
{
  as.double(pROC::auc(roc(y,y_hat,quiet = T)))
}


#### Main function for fitting logistic regression model ####
logistic.CD <- function(x,y,y2=NULL,penalty="lasso",
                        lambda=NULL,gamma=NULL,
                        beta_support=NULL, standardize=T, maxit=100,tol=1e-4,
                        x_val=NULL,y_val=NULL,beta.oracle=NULL,
                        x_test=NULL,y_test=NULL,
                        best=F,metric="bic")
{
  #### function used to fit a penalized logistic regression using coordinate descent ####
  # Args:
  #  x: covariates matrix (n, p) 
  #  y: resposne vector (n, ) used for model training
  #  y2: resposne vector (n, ) used for model evaluation
  #      In most of cases, y2=y.
  #  penalty: one of "lasso","MCP","SCAD"
  #  lambda: can be a sequence of tuning parameters
  #  gamma: gamma parameter for MCP or SCAD 
  #  beta_support: indicator vector (p,), indicating variables must stay in model (so no penalize is imposed in them)
  #  standardize: logic, whether to standarize covairates before model training
  #  maxit: maximum number of iterations for coordinate descent
  #  tol: tolerance for convergence
  #  x_val, y_val: validation set 
  #  x_test, y_test: test set
  #  beta.oracle: true coefficients vector  
  #  best: whether to output the best model performance (regarding different tunings) only, otherwise, results for all tunings are output
  #  metric: metric to select the best model, one of "bic","aic","loglike_val","AUC_val"

  # Return:
  #  result: a list containing the model information and model evaluation metrics (variable selection, estimation, prediction performance etc.)
  
  #### Initialization ####
  n = dim(x)[1]
  p = dim(x)[2]
  
  # standardize x
  if(standardize){
    centerx = apply(x,2,mean)
    stdx = apply(x,2,sd)
    x = (x - outer(rep(1,n),centerx,"*"))/outer(rep(1,n),stdx)
  }
  
  if(is.null(beta_support)){
    beta_support = rep(F,p)
  }
  
  if(is.null(y2)){
    y2 = y
  }
  
  if(is.null(gamma)){
    if(penalty=="SCAD"){
      gamma = 3.7
    }else if(penalty=="MCP"){
      gamma = 3
    }
  }
  
  if(is.null(lambda)){
    lambda = seq(1,0.01,length.out = 20)
  }
  lambda = sort(lambda,decreasing = T)
  n.lambda = length(lambda)
  
  beta = matrix(NA,p+1,n.lambda) # matrix containing result(\beta_0 and \beta) for each lambda
  convergence = rep(T,n.lambda)
  iterations = rep(NA,n.lambda)
  
  #### Iteration ####
  for(i in 1:n.lambda){
    lambda.tmp = lambda[i]
    if((i>1)&(convergence[max(i-1,1)])){ 
      beta0 = as.vector(beta[,i-1])  
      #beta0 = c(-3,runif(p,-0.5,0.5))
    }else{
      beta0 = c(-3,runif(p,-0.5,0.5))
    }
    for(iter in 1:maxit){
      beta_pre = beta0 # beta_pre is used to calculate stoping criterion latter
      ## update for \beta
      j_order_seq = order(beta_support,decreasing = T)
      for(j_order in j_order_seq){
        j = (1:p)[j_order]
        P = expit(x%*%beta0[-1] + beta0[1])
        w_j = sum((y-P)*x[,j])/n
        #v_j = sum((P*(1-P)*x[,j]^2))/n
        p_factor = max(mean(P*(1-P)),0.01) # for stable of the algorithm
        v_j = p_factor*mean(x[,j]^2) 
        z_j = beta0[j+1]+w_j/v_j
        ## thresholding according to "Adaptive rescaling" in paper
        ## "COORDINATE DESCENT ALGORITHMS FOR NONCONVEX PENALIZED REGRESSION, WITH APPLICATIONS TO BIOLOGICAL FEATURE SELECTION"
        if(beta_support[j]){
          beta_j_up = z_j
        }else{
          if(penalty=="lasso"){
            beta_j_up = soft.thresholding(z_j,lambda.tmp/v_j)
          }else if(penalty=="SCAD"){
            #print(sprintf("%.5f,%.5f,%.5f,%d,%.5f",z_j,v_j,lambda[i],i,lambda.tmp))
            beta_j_up = scad.thresholding(z_j*v_j,lambda.tmp,gamma) / v_j
            #print(sprintf("%.5f,%.5f",z_j,v_j))
          }else if(penalty=="MCP"){
            beta_j_up = mcp.thresholding(z_j*v_j,lambda.tmp,gamma) / v_j
          }
        }
        beta0[j+1] = beta_j_up
      }
      
      ## update for \beta_0
      P = expit(x%*%beta0[-1] + beta0[1])
      p_factor = max(mean(P*(1-P)),0.01)
      #P = sapply(expit(x%*%beta0[-1] + beta0[1]),function(x) max(0.001,min(x,0.999)))
      beta0[1] = beta0[1] + mean(y-P)/p_factor
      beta0[1] = max(min(beta0[1],1e3),-1e3)
      if(sum((beta0[-1]-beta_pre[-1])^2)/p<tol){
        break
      }
    }
    beta[,i] = beta0
    iterations[i] = iter
    if(iter==maxit){
      convergence[i] = F
      #cat("May not convergence: lambda= ",lambda.tmp,"\n")
    }
  }
  
  #### post processing ####
  if(standardize){
    beta_out = beta
    for(i in 1:n.lambda){
      beta_out[-1,i] = beta[-1,i] / stdx
      beta_out[1,i] = beta[1,i] - sum(beta[-1,i]*centerx/stdx)
    }
  }else{
    beta_out = beta
  }
  
  #### output ####
  ACURACCY_val = rep(NA,n.lambda)
  AUC_val = rep(NA,n.lambda)
  loglike_val = rep(NA,n.lambda)
  AUC_test = rep(NA,n.lambda)
  loglike_test = rep(NA,n.lambda)
  loglike = rep(NA,n.lambda)
  aic = rep(NA,n.lambda)
  bic = rep(NA,n.lambda)
  FNR = rep(NA,n.lambda)
  FPR = rep(NA,n.lambda)
  FR = rep(NA,n.lambda)
  CNZ = rep(NA,n.lambda)
  INZ = rep(NA,n.lambda)
  MCC = rep(NA,n.lambda)
  FDR = rep(NA,n.lambda)
  BETA_MSE = rep(NA,n.lambda)
  BETA_BIAS = rep(NA,n.lambda) #L1 norm bias of an estimator for the nonzero coefficients;
  
  for(i in 1:n.lambda){
    beta.tmp = as.vector(beta[,i])
    k = sum(beta.tmp[-1]!=0)
    y_hat = as.vector(expit(x%*%beta.tmp[-1]+beta.tmp[1]))
    loglike[i] = loglike.logistic(y2,y_hat)
    bic[i] = -2*n*loglike[i]+k*log(n)
    aic[i] = -2*n*loglike[i]+k*2
    
    if(!is.null(beta.oracle)){
      info = cal.coeff_metrics(beta.oracle[-1],beta.tmp[-1])
      FNR[i] = info$FNR
      FPR[i] = info$FPR
      FR[i] = info$FR
      CNZ[i] = info$CNZ
      INZ[i] = info$INZ
      MCC[i] = info$MCC
      FDR[i] = info$FDR
      BETA_MSE[i] = info$BETA_MSE
      BETA_BIAS[i] = info$BETA_BIAS
    }
    if(!is.null(x_val)){
      y_val_hat = as.vector(expit(x_val%*%beta.tmp[-1]+beta.tmp[1]))
      ACURACCY_val[i] = sum((y_val_hat>=0.5)==(y_val>=0.5))/length(y_val)
      AUC_val[i] = cal.auc(y_val,y_val_hat)
      loglike_val[i] = loglike.logistic(y_val,y_val_hat)
    }
    if(!is.null(x_test)){
      y_test_hat = as.vector(expit(x_test%*%beta.tmp[-1]+beta.tmp[1]))
      AUC_test[i] = cal.auc(y_test,y_test_hat)
      loglike_test[i] = loglike.logistic(y_test,y_test_hat)
    }
  }
  
  if(metric=="bic"){
    metric_val = bic+(1-convergence)*1e10
  }else if(metric=="aic"){
    metric_val = aic+(1-convergence)*1e10
  }else if(metric=="loglike_val"){
    metric_val = -loglike_val+(1-convergence)*1e10
  }else if(metric=="AUC_val"){
    metric_val = -AUC_val+(1-convergence)*1e10
  }
  
  if(best){
    best_idx = which.min(metric_val)
    return(list(beta=beta_out[,best_idx],lambda=lambda[best_idx],
                loglike = loglike[best_idx], aic=aic[best_idx],bic=bic[best_idx],
                FNR = FNR[best_idx],FPR=FPR[best_idx],FR=FR[best_idx],CNZ=CNZ[best_idx],INZ=INZ[best_idx],MCC=MCC[best_idx],FDR=FDR[best_idx],
                BETA_MSE=BETA_MSE[best_idx],BETA_BIAS=BETA_BIAS[best_idx],
                loglike_val=loglike_val[best_idx],AUC_val=AUC_val[best_idx],ACURACCY_val=ACURACCY_val[best_idx],
                loglike_test=loglike_test[best_idx],AUC_test=AUC_test[best_idx],
                metric_val=metric_val[best_idx]
      )
    )
    
  }
  return(list(beta=beta_out,lambda=lambda,
              convergence=convergence,iterations=iterations,
              loglike = loglike, aic=aic,bic=bic,
              FNR = FNR,FPR=FPR,FR=FR,CNZ=CNZ,INZ=INZ,MCC=MCC,FDR=FDR,
              BETA_MSE=BETA_MSE,BETA_BIAS=BETA_BIAS,
              loglike_val=loglike_val,AUC_val=AUC_val,ACURACCY_val=ACURACCY_val,
              loglike_test=loglike_test,AUC_test=AUC_test,
              metric_val=metric_val
              )
         )
}


#### CV ####
cv.logistic <- function(x,y,y2=NULL, penalty="lasso",
                        lambda=NULL, gamma=NULL,
                        beta_support=NULL, standardize=T, maxit=100,tol=1e-5,
                        nfolds=4,fold=NULL,
                        metric = "loglike_val",onestd = F,
                        x_val=NULL,y_val=NULL,beta.oracle=NULL
                        )
{
  #### Cross validation to find the best tuning for logistic regression model, ####
  #### and used the best tuning to fit a model based on full dataset
  
  # Args:
  #  x: covariates matrix (n, p) 
  #  y: resposne vector (n, ) used for model training
  #  y2: resposne vector (n, ) used for model evaluation
  #      In most of cases, y2=y.
  #  penalty: one of "lasso","MCP","SCAD"
  #  lambda: can be a sequence of tuning parameters
  #  gamma: gamma parameter for MCP or SCAD 
  #  beta_support: indicator vector (p,), indicating variables must stay in model (so no penalize is imposed in them)
  #  standardize: logic, whether to standarize covairates before model training
  #  maxit: maximum number of iterations for coordinate descent
  #  tol: tolerance for convergence
  
  #  nfolds: number of folds for CV
  #  fold: a vector (n,) indicating the fold each sample belongs to (can be NULL, then a new n-fold split is created)
  #  metric: metric for CV to select the best model, one of "loglike_val","AUC_val"
  #  onestd: logic, whether "mean + 1 std" strategy is used to select the best model (default False) 
  
  #  x_val, y_val: validation set 
  #  beta.oracle: true coefficients vector  
  
  # Return:
  #  result: a list containing the model information and model evaluation metrics (variable selection, estimation, prediction performance etc.)
  n = dim(x)[1]
  if(is.null(y2)){
    y2 = y
  }
  if(is.null(fold)){
    idx0 = which(y<0.5)
    idx1 = which(y>0.5)
    tmp0 = createFolds(1:length(idx0),k=nfolds,list=TRUE,returnTrain=F)
    tmp1 = createFolds(1:length(idx1),k=nfolds,list=TRUE,returnTrain=F)
    fold = rep(NA,n)
    for(i in 1:nfolds){
      idx = c(idx0[tmp0[[i]]],idx1[tmp1[[i]]])
      fold[idx] = i
    }
  }else{
    nfolds = length(unique(fold))
  }
  
  if(is.null(lambda)){
    lambda = seq(1,0.01,length.out = 20)
  }
  n.lambda = length(lambda)
  
  metric_cv = matrix(NA,nfolds,n.lambda)
  loglike_val_cv = matrix(NA,nfolds,n.lambda)
  AUC_val_cv = matrix(NA,nfolds,n.lambda)
  for(i in 1:nfolds){
    idx_val = which(fold==i)
    x_train = x[-idx_val,]
    y_train = y[-idx_val]
    x_val_cv = x[idx_val,]
    y_val_cv = y2[idx_val]
    m = logistic.CD(x_train,y_train,maxit=maxit,tol=tol,penalty=penalty,
                    beta_support=beta_support,
                    lambda=lambda,gamma=gamma,standardize=standardize,
                    x_val=x_val_cv,y_val=y_val_cv,beta.oracle=beta.oracle,metric=metric,
                    )
    metric_cv[i,] = m$metric_val
    
    loglike_val_cv[i,] = m$loglike_val
    AUC_val_cv[i,] = m$AUC_val
  }
  
  metric_cv_mean = apply(metric_cv,2,mean)
  if(onestd){
    metric_cv_used = metric_cv_mean + metric_val_cv_std # min is better
  }else{
    metric_cv_used = metric_cv_mean
  }
  best_idx = which.min(metric_cv_used)
  
  #### use best lambda to fit the whole model ####
  lambda_best = lambda[best_idx]
  model = logistic.CD(x,y,maxit=maxit,tol=tol,penalty=penalty,
                      beta_support=beta_support,
                      lambda=lambda_best,gamma=gamma,standardize=standardize,
                      x_val=x_val,y_val=y_val,beta.oracle=beta.oracle,metric=metric
  )
  
  return(list(lambda=lambda,
              metric_cv_used=metric_cv_used,
              metric_cv_used_best = metric_cv_used[best_idx],
              loglike_val_cv=apply(loglike_val_cv,2,mean)[best_idx],
              AUC_val_cv=apply(AUC_val_cv,2,mean)[best_idx],
              best=model))
}

  
#### Main function for fitting the proposed multi-prior transfer learning model ####
transfer.logistic <- function(x, y, y_hat=NULL, eta=NULL,penalty="lasso", 
                              lambda=NULL,gamma=NULL,
                              standardize=T, maxit=100, tol=1e-5,
                              x_val=NULL,y_val=NULL,x_test=NULL,y_test=NULL,beta.oracle=NULL,
                              best=F,metric="bic")
{
  #### function used to fit the proposed model using coordinate descent ####
  # Args:
  #  x: covariates matrix (n, p) 
  #  y: resposne vector (n, ) used for model training
  #  y_hat: matrix for K prior y (p, K)
  #  eta: weight vector (K,) for each prior (here, the eta is set to be (global weight)*(relative weights) in paper)
  
  #  penalty: one of "lasso","MCP","SCAD"
  #  lambda: can be a sequence of tuning parameters
  #  gamma: gamma parameter for MCP or SCAD 
  #  beta_support: indicator vector (p,), indicating variables must stay in model (so no penalize is imposed in them)
  #  standardize: logic, whether to standarize covairates before model training
  #  maxit: maximum number of iterations for coordinate descent
  #  tol: tolerance for convergence
  #  x_val, y_val: validation set 
  #  x_test, y_test: test set
  #  beta.oracle: true coefficients vector  
  #  best: whether to output the best model performance (regarding different tunings) only, otherwise, results for all tunings are output
  #  metric: metric to select the best model, one of "bic","aic","loglike_val","AUC_val"
  
  # Return:
  #  result: a list containing the model information and model evaluation metrics (variable selection, estimation, prediction performance etc.)
  
  
  #### compute surrogate y (y_hat) ####
  eta = eta
  y2 = y
  if(!is.null(y_hat)){
    y = (y_hat%*%eta + y) / (1+sum(eta)) #surrogate y
  }
  
  #### model fitting ####
  # fit the entire path 
  model = logistic.CD(x,y,y2=y2,maxit=maxit,tol=tol,penalty=penalty,
                      lambda=lambda,gamma=gamma,standardize=standardize,
                      x_val=x_val,y_val=y_val,x_test=x_test,y_test=y_test,beta.oracle=beta.oracle,
                      best=best,metric = metric)
  
  #### Output ####
  model[["eta"]] = eta
  return(model)
}

#### CV ####
cv.transfer.logistic <- function(x, y, y_hat=NULL, v=1, eta0=NULL,penalty="lasso",
                                 lambda=NULL, gamma=NULL, 
                                 standardize=T, maxit=100, tol=1e-5,
                                 nfolds=4,fold=NULL,
                                 metric = "loglike_val",onestd = F,
                                 x_val=NULL,y_val=NULL,beta.oracle=NULL
)
{
  #### Cross validation to find the best tuning for the proposed method, ####
  #### and used the best tuning to fit a model based on full dataset
  
  # Args:
  #  x: covariates matrix (n, p) 
  #  y: resposne vector (n, ) used for model training
  #  y_hat: matrix for K prior y (p, K)
  #  v: prior global weight
  #  eta0: weight vector (K,) for each prior (relative weights)
  
  #  penalty: one of "lasso","MCP","SCAD"
  #  lambda: can be a sequence of tuning parameters
  #  gamma: gamma parameter for MCP or SCAD 
  #  beta_support: indicator vector (p,), indicating variables must stay in model (so no penalize is imposed in them)
  #  standardize: logic, whether to standarize covairates before model training
  #  maxit: maximum number of iterations for coordinate descent
  #  tol: tolerance for convergence
  
  #  nfolds: number of folds for CV
  #  fold: a vector (n,) indicating the fold each sample belongs to (can be NULL, then a new n-fold split is created)
  #  metric: metric for CV to select the best model, one of "loglike_val","AUC_val"
  #  onestd: logic, whether "mean + 1 std" strategy is used to select the best model (default False) 
  
  #  x_val, y_val: validation set 
  #  beta.oracle: true coefficients vector  
  
  # Return:
  #  result: a list containing the model information and model evaluation metrics (variable selection, estimation, prediction performance etc.)
  
  y2 = y
  n_v = length(v)
  best_metric = Inf
  result = NULL
  for(i in 1:n_v){
    eta = eta0*v[i]
    if(!is.null(y_hat)){
      y = (y_hat%*%eta + y) / (1+sum(eta)) #surrogate y
    }
    result_tmp = cv.logistic(x,y,y2=y2,penalty=penalty,lambda=lambda,nfolds=nfolds,fold=fold,
                             metric = metric,onestd = onestd,
                             maxit=maxit,tol=tol,gamma=gamma,standardize=standardize,
                             x_val=x_val,y_val=y_val,beta.oracle=beta.oracle)
    result_tmp$best[["eta"]] = eta
    if(result_tmp$metric_cv_used_best<best_metric){
      result = result_tmp
      best_metric = result_tmp$metric_cv_used_best
    }
  }
  return(result)
}

#### calculation of weights ####
weight.prior <- function(loglike_val, v=1)
{
  # Args:
  #  loglike_val: validation log-likelihood for each prior (K,)  
  #  v: glocal prior weight
  t = sum(exp(loglike_val))
  eta = v * exp(loglike_val) / t
  return(eta)
}




#### cacluate the beta_support from prior datasets (used in real data analysis only) ####
gen_beta_support <- function(X,Y,
                             lambda.seq=round(exp(seq(log(0.3), log(0.002),length.out=8)),digits=5))
{
  model = logistic.CD(X,Y,penalty="SCAD", standardize = F,lambda = lambda.seq,
                      best=T,metric = "bic"
  )
  beta_support = (model$beta[-1]!=0)
  return(list(beta_support=beta_support,model=model))
}