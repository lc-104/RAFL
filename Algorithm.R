#####This file contains the main functions we used to implement the RAFL#####
#####including functions for Algorithm 1, 2, and 3#####


#####Robust loss functions used in the paper#####

#Tukey biweight loss function
tukey_loss=function(x, M){
  y = x^2/2-x^4/(2*M^2)+x^6/(6*M^4)
  y = y*(abs(x)<M)+M^2/6*(abs(x)>=M)
  return(y)
}

#first derivative of the Tukey biweight loss function
tukey_psi=function(x, M){
  y = x*(1-(x/M)^2)^2;
  y = y*(abs(x)<M)+0*(abs(x)>=M)
  return(y)
}

#Huber's loss function
huber_loss=function(x, M){
  y = x^2/2
  y = y*(abs(x)<M)+(M*(abs(x)-0.5*M))*(abs(x)>=M)
  return(y)
}

#first derivative of the Huber's loss function
huber_psi=function(x, M){
  y = x
  y = y*(abs(x)<M)+M*sign(x)*(abs(x)>=M)
  return(y)
}

#####apg algorithm to solve the standard lasso
###lambda1: lasso penalty =0 for ridge only
###lambda2: ridge penalty =0 for lasso only


apglasso=function(X, Y, type, M, lambda1, lambda2, weight=1, sig=1, beta0=NULL, MAX_ITERS=50, EPS = 1e-4, QUIET=FALSE){
  
  n=nrow(X)
  p=ncol(X)
  Ymean=mean(Y)
  Xmean=colMeans(X)
  Y=scale(Y,scale=FALSE)
  X=scale(X,scale=FALSE)
 
  #gradient and proximal operator
  
  if(type=="square"){
    #square loss
    #f <- function(beta){sum((X%*%beta-Y)^2)/n}
    gradf <- function(beta){2*t(X)%*%(X%*%beta-Y)/n+2*lambda2*beta}
  } else if(type=="huber"){
    #huber loss
    #f <- function(beta){2*sum(huber_loss((X%*%beta-Y)/sig,M))/n}
    gradf <- function(beta){2*t(X)%*%huber_psi((X%*%beta-Y)/sig,M)/sig/n+2*lambda2*beta}
  } else if(type=="tukey"){
    #tukey loss
    #f <- function(beta){2*sum(tukey_loss((X%*%beta-Y)/sig,M))/n}
    gradf <- function(beta){2*t(X)%*%tukey_psi((X%*%beta-Y)/sig,M)/sig/n+2*lambda2*beta}
  } else{stop("loss function not defined")
  }
  
  ###non-smooth part, lasso penalty, proximal soft thresholding
  #g <- function(beta) {sum(abs(beta)*(lambda*weight))}
  proxg <- function(beta, t) {sign(beta)*(sapply(abs(beta) - t*(lambda1*weight),
                                                 FUN=function(x) {max(x,0)})) }
  
  ###initial values
  if(is.null(beta0)){
    beta0 <- rep(0,p)
  }
  
  sol=apg(gradf, proxg, p, opts=list(X_INIT=beta0, MAX_ITERS=MAX_ITERS, EPS = EPS, QUIET=QUIET))
  
  mu=as.vector(Ymean-t(sol$x)%*%Xmean)
  
  scale=mad(Y-X%*%sol$x-mu)
  
  return(list(beta=sol$x, mu=mu, scale=scale))
}



#####Algorithm 1: APG for sub-problem of the beta update#####

apgbeta=function(X, Y, type=c("square","huber","tukey"), M, sig, tau, D, theta, u, lambda, weight=1, X_INIT, MAX_ITERS=50, EPS = 1e-4, QUIET=FALSE){
  
  
  n=nrow(X)
  p=ncol(X)
  
  #objective functions with extra admm term
  #admm term  
  #h <- function(beta){tau/2*norm_vec(D%*%beta-theta+u)}  
  #gradh <- function(beta){tau*t(D)%*%(D%*%beta-theta+u)}
  
  if(type=="square"){
  #square loss
  #f <- function(beta){sum((X%*%beta-Y)^2)/n}
  gradf <- function(beta){2*t(X)%*%(X%*%beta-Y)/n+tau*t(D)%*%(D%*%beta-theta+u)}
  } else if(type=="huber"){
  #huber loss
  #f <- function(beta){2*sum(huber_loss((X%*%beta-Y)/sig,M))/n}
  gradf <- function(beta){2*t(X)%*%huber_psi((X%*%beta-Y)/sig,M)/sig/n+tau*t(D)%*%(D%*%beta-theta+u)}
  } else if(type=="tukey"){
  #tukey loss
  #f <- function(beta){2*sum(tukey_loss((X%*%beta-Y)/sig,M))/n}
  gradf <- function(beta){2*t(X)%*%tukey_psi((X%*%beta-Y)/sig,M)/sig/n+tau*t(D)%*%(D%*%beta-theta+u)}
  } else{stop("loss function not defined")
  }
  
  ###non-smooth part, lasso penalty, proximal soft thresholding
  #g <- function(beta) {sum(abs(beta)*(lambda*weight))}
  proxg <- function(beta, t) {sign(beta)*(sapply(abs(beta) - t*(lambda*weight),
                                                    FUN=function(x) {max(x,0)})) }
  
  
  
  sol=apg(gradf, proxg, p, opts=list(X_INIT=X_INIT, MAX_ITERS=MAX_ITERS, EPS = EPS, QUIET=QUIET))
  
  return(list(beta =sol$x))
}


#####Algorithm 2: APG for sub-problem the theta update#####

apgtheta=function(tau, D, beta, u, lambda, X_INIT, MAX_ITERS=50, EPS = 1e-4, QUIET=FALSE){
  

  #objective functions with extra admm term
  #admm term  
  #h <- function(theta){tau/2*norm_vec(D%*%beta-theta+u)}  
  gradh <- function(theta){tau*(theta-D%*%beta-u)}
  
  ###non-smooth part, lasso penalty, proximal soft thresholding
  #g <- function(theta) {sum(abs(theta)*(lambda))}
  proxg <- function(theta, t) {sign(theta)*(sapply(abs(theta) - t*(lambda),
                                                 FUN=function(x) {max(x,0)})) }
  
  sol=apg(gradh, proxg, length(X_INIT), opts=list(X_INIT=X_INIT, MAX_ITERS=MAX_ITERS, EPS = EPS, QUIET=QUIET))
  
  return(list(theta =sol$x))
}




#####Algorithm 3: ADMM embedded with two APG algorithm to solve the robust adaptive fused lasso #####
###input:
#X: design matrix
#Y: response matrix
#lambda1 fused lasso penalty
#lambda2 sparsity penalty
#type: type of the loss function

raflasso=function(X, Y, type=c("square","huber","tukey"), M, sig, D, lambda1, lambda2, weight=1, scale_admm=1, beta0=NULL, MAX_ITERS=50, EPS = 1e-4, QUIET=TRUE, max_iters_admm=50, convergence_admm = 1e-4, change_scale=FALSE){
  
  n=nrow(X)
  p=ncol(X)
  Ymean=mean(Y)
  Xmean=colMeans(X)
  Y=scale(Y,scale=FALSE)
  X=scale(X,scale=FALSE)
  
  #scale_admm=1
  
  ###initial values
  if(is.null(beta0)){
    beta0 <- rep(0,p)
  }
  
  beta=beta0
  theta=D%*%beta
  u=rep(0,dim(D)[1])
  
  dualNorm=NULL
  primalNorm=NULL
  
  for(i in 1:max_iters_admm){
    ##update beta
    betaupdate=apgbeta(X, Y, type=type, M=M, sig=sig, tau=scale_admm, D=D, theta=theta, u=u, lambda=lambda2, weight=weight, X_INIT=beta, MAX_ITERS=MAX_ITERS, EPS = EPS, QUIET=QUIET)
      
    beta=betaupdate$beta
    thetaOld = theta
    
    
    ##update theta
    thetaupdate= apgtheta(tau=scale_admm, D=D, beta=beta, u=u, lambda=lambda1, X_INIT=thetaOld, MAX_ITERS=MAX_ITERS, EPS = EPS, QUIET=QUIET)
    theta=thetaupdate$theta
    
    ##update scaled dual variables u
    dualResNorm = norm_vec(theta-thetaOld)/ max(1,norm_vec(thetaOld))
    primalRes = D%*%beta - theta
    primalResNorm = norm_vec(primalRes)/ max(1,norm_vec(theta))
    u = u + primalRes
    
    
    dualNorm=c(dualNorm,dualResNorm)
    primalNorm=c(primalNorm,primalResNorm)
    
    
    if (primalResNorm < convergence_admm
        & dualResNorm < convergence_admm) {
      break
    }
    
    # update ADMM scale parameter u
    if (change_scale==TRUE){
      if (primalResNorm/dualResNorm>10){
        scale_admm=scale_admm/2
        u = u/2     
      } else if (primalResNorm/dualResNorm<0.1) {
        scale_admm=scale_admm*2
        u = 2*u 
      }
    }
  }

  mu=Ymean-t(beta)%*%Xmean
  
  return(list(beta =beta, theta=theta, mu=as.vector(mu), nit=i, dualNorm=dualNorm, primalNorm=primalNorm, scale_admm=scale_admm))
}




#####fucntion to select tuning parameters using TMSPE#####

raflasso.cv=function(X, Y, K, trim=0.8, lambda1.min=0.01, lambda2.min=0.01, lambda1.max=1, lambda2.max=1, nlambda1=5, nlambda2=5,  type=c("square","huber","tukey"), M, sig, D, weight=1, 
                     scale_admm=1, beta0=NULL, MAX_ITERS=50, EPS = 1e-6, QUIET=TRUE, max_iters_admm=50, convergence_admm = 1e-4, change_scale=FALSE){
  n=nrow(X)
  p=ncol(X)
  Data=cbind(Y,X)
  #Randomly shuffle the data
  sData<-Data[sample(n),]
  
  #Create K equally size folds
  folds <- cut(seq(1,nrow(sData)),breaks=K,labels=FALSE)
  testData=list()
  trainData=list()
  
  #Perform K fold cross validation
  for(i in 1:K){
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData[[i]] <- sData[testIndexes, ]
    trainData[[i]] <- sData[-testIndexes, ]
  }
  
  lambda1seq=seq(lambda1.min,lambda1.max,length.out=nlambda1)
  lambda2seq=seq(lambda2.min,lambda2.max,length.out=nlambda2)
  lambda=expand.grid(lambda1=lambda1seq,lambda2=lambda2seq)
  
  criteria=vector()

  for (k in 1:nrow(lambda)){
    lambda1=lambda[k,1]
    lambda2=lambda[k,2]
    
    cv.error=NULL
    for (i in 1:K){
      train=trainData[[i]]
      test=testData[[i]]
      sol=raflasso(train[,-1], train[,1], type=type, M=M, sig=sig, D=D, weight=weight, lambda1=lambda1, lambda2=lambda2,
                   scale_admm=scale_admm, beta0=beta0, MAX_ITERS=MAX_ITERS, EPS = EPS, QUIET=QUIET, max_iters_admm=max_iters_admm, convergence_admm = convergence_admm , change_scale= change_scale)
      
      ###calculate prediction error
      error=test[,1]-test[,-1]%*%sol$beta-as.vector(sol$mu)
      errorsq=sort(error^2)[1:as.integer(trim*(length(error^2)))]
      cv.error=c(cv.error, sum(errorsq))
    }  
    criteria[k]=mean(cv.error)
  }
  mink=which.min(criteria)
  lambda1=lambda[mink,1]
  lambda2=lambda[mink,2]
  
  return(list(error.CV =criteria, lambda1=lambda1, lambda2=lambda2))
}


#####main function to perform APG algorithm##### modified from https://github.com/jpvert/apg/blob/master/R/apg.R

apg <- function(grad_f, prox_h, dim_x, opts) {
  
  # Set default parameters
  X_INIT <- numeric(dim_x) # initial starting point
  USE_RESTART <- TRUE # use adaptive restart scheme
  MAX_ITERS <- 100 # maximum iterations before termination
  EPS <- 1e-6 # tolerance for termination
  ALPHA <- 1.01 # step-size growth factor
  BETA <- 0.5 # step-size shrinkage factor
  QUIET <- FALSE # if false writes out information every 100 iters
  GEN_PLOTS <- TRUE # if true generates plots of norm of proximal gradient
  USE_GRA <- FALSE # if true uses UN-accelerated proximal gradient descent (typically slower)
  STEP_SIZE = NULL # starting step-size estimate, if not set then apg makes initial guess
  FIXED_STEP_SIZE <- FALSE # don't change step-size (forward or back tracking), uses initial step-size throughout, only useful if good STEP_SIZE set
  
  # Replace the default parameters by the ones provided in opts if any
  for (u in c("X_INIT","USE_RESTART", "MAX_ITERS", "EPS", "ALPHA", "BETA", "QUIET", "GEN_PLOTS", "USE_GRA", "STEP_SIZE", "FIXED_STEP_SIZE")) {
    eval(parse(text=paste('if (exists("',u,'", where=opts)) ',u,' <- opts[["',u,'"]]',sep='')))
  }
  
  # Initialization
  x <- X_INIT
  y <- x
  g <- grad_f(y)
  if (norm_vec(g) < EPS) {
    return(list(x=x,t=0))
  }
  theta <- 1
  
  # Initial step size
  if (is.null(STEP_SIZE)) {
    
    # Barzilai-Borwein step-size initialization:
    t <- 1 / norm_vec(g)
    x_hat <- x - t*g
    g_hat <- grad_f(x_hat)
    t <- abs(sum( (x - x_hat)*(g - g_hat)) / (sum((g - g_hat)^2)))
  } else {
    t <- STEP_SIZE
  }
  
  # Main loop
  for (k in seq(MAX_ITERS)) {
    
    if (!QUIET && (k %% 100==0)) {
      message(paste('iter num ',k,', norm(tGk): ',err1,', step-size: ',t,sep=""))
    }
    
    x_old <- x
    y_old <- y
    
    # The proximal gradient step (update x)
    x <- prox_h( y - t*g, t)
    
    # The error for the stopping criterion
    err1 <- norm_vec(y-x) / max(1,norm_vec(x))
    if (err1 < EPS) break
    
    # Update theta for acceleration
    if(!USE_GRA)
      theta <- 2/(1 + sqrt(1+4/(theta^2)))
    else
      theta <- 1
    end
    
    # Update y
    if (USE_RESTART && sum((y-x)*(x-x_old))>0) {
      x <- x_old
      y <- x
      theta <- 1
    } else {
      y <- x + (1-theta)*(x-x_old)
    }
    
    # New gradient
    g_old <- g
    g <- grad_f(y)
    
    # Update stepsize by TFOCS-style backtracking
    if (!FIXED_STEP_SIZE) {
      t_hat <- 0.5*sum((y-y_old)^2)/abs(sum((y - y_old)*(g_old - g)))
      t <- min( ALPHA*t, max( BETA*t, t_hat ))
    }
  }
  if (!QUIET) {
    message(paste('iter num ',k,', norm(tGk): ',err1,', step-size: ',t,sep=""))
    if (k==MAX_ITERS) message(paste('Warning: maximum number of iterations reached'))
    message('Terminated')
  }
  
  # Return solution and step size
  return(list(x=x,t=t))
}


norm_vec <- function(x) {
  sqrt(sum(x^2))
}

MMRidge<- function (X,Y,lambda){
  sol=apglasso(X, Y, type="huber", 5, lambda1=0, lambda2=lambda, weight=1, sig=1, beta0=NULL, MAX_ITERS=50, EPS = 1e-4, QUIET=FALSE)
  return(sol)  
}
