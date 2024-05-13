library(MASS)

data_generator1 <- function(n=200,p=100,n_val=200,n_test=200,
                            beta=c(0,
                                   -0.5,-2,0.5,2,-1.5,
                                   1,2,-1.5,2,-2,
                                   1,1.5,-2,1,1.5,rep(0,p-15)),
                            discret = NULL,
                            rou = 0.5, seed=NULL)
{
  #### data generating function ####
  # Args:
  #  n: sample size
  #  p: num of covariates
  #  n_val: num of validation samples
  #  n_test: num of test samples
  #  beta: true beta coefficients
  #  discret: the idxs set of discret variables
  #  rou: parameter for the AR structure predictor covariance matrix 
  #  seed: random seeed
  
  # Return:
  #  result: a list containing the generated datasets
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  if(is.null(beta)){
    beta = c(0,runif(15,-3,3),rep(0,p-15))
  }
  if(is.null(discret)){
    # set 40% of variables to be binary (first 2 variables of every 5 variabels)
    discret = c(seq(1,p,5),seq(2,p,5))
  }
  mu = rep(0,p)
  Sigma = rou^abs(outer(1:p,1:p,"-"))
  x_all = mvrnorm(n+n_val+n_test,mu,Sigma)
  x_all[,discret] = sign(x_all[,discret]>0)
  prob_all = 1/(1+exp(-(beta[1]+x_all%*%beta[-1])))
  y_all = rbinom(n+n_val+n_test,1,prob_all)
  
  x = x_all[1:n,]
  y = y_all[1:n]
  prob = prob_all[1:n]
  
  x_val = x_all[(n+1):(n+n_val),]
  y_val = y_all[(n+1):(n+n_val)]
  
  x_test = x_all[(n+n_val+1):(n+n_val+n_test),]
  y_test = y_all[(n+n_val+1):(n+n_val+n_test)]
  
  return(list(x=x,y=y,x_val=x_val,y_val=y_val,x_test=x_test,
              y_test=y_test,beta=beta,prob=prob)
         )
}







