#####This file contains the codes used to run the simulation study Example 1 under all three scenarios.#####
library(mvtnorm)
library(MASS)
library(Matrix)
library(genlasso)
library(glmnet)
library(EvaluationMeasures)
library(pdfCluster)
source('Algorithm.R')


###level of contamination
###setting for scenario 1: lv.X=0,lv.Y=0
###setting for scenario 2: lv.X=0,lv.Y=0.1
###setting for scenario 3: lv.X=0.1,lv.Y=0.1

lv.X=0.1
lv.Y=0.1

###set up parameters
n=100
p=50
var=1
s=20
r=1
h=0.1

###run M simulations
M=100
normball=matrix(NA,M,7)
MSPEall=matrix(NA,M,7)
TPRall=matrix(NA,M,7)
TNRall=matrix(NA,M,7)
randall=matrix(NA,M,7)

###generate design matrix
set.seed(1)

corx=0
covX=diag(p)
for(i in 1:p) for(j in 1:p)
  covX[i,j]=corx^abs(i-j)

X=rmvnorm(n,mean=rep(0,p),sigma = covX)
X.ori=X

###add contamination on X
cont=matrix(0,n,p)
cont.id=cbind(rep(sample(n,ceiling(n*lv.X)),ceiling(p/50)),sample(p,ceiling(n*lv.X)*ceiling(p/50),replace=TRUE))
cont[cont.id]=rnorm(nrow(cont.id), 10, 1)
X=X+cont

X.test=rmvnorm(n,mean=rep(0,p),sigma = covX)


for (k in 1:M){

###generate data
set.seed(k)

###generate B 
coef=c(rep(c(-4*r,-2*r,2*r,4*r), each=s/4),rep(0, p-s))
group=c(rep(1:4,each=s/4),rep(5,p-s))
b0=coef
trueb=as.numeric(as.vector(b0)!=0)

E=rnorm(n,0,var)

Y.ori=X.ori%*%b0 + E

###add contamination on errors
cont=rep(0,n)
cont.id=sample(n,ceiling(n*lv.Y))
cont[cont.id]=rnorm(length(cont.id), 10, 1)

E.cont=E+cont
Y=X.ori%*%b0 + E.cont

###testing set
E.test=rnorm(n,0,var)
Y.test=X.test%*%b0+E.test


####lasso and adaptive lasso
reg.fit<-lm(Y~X)
weights=1/abs(reg.fit$coef[-1])^0.5

cvfit=cv.glmnet(x=X,y=Y,alpha=1,nfolds=5)

pred0=predict(cvfit, newx=X.test, s="lambda.min")
norm0=sum(coef(cvfit,s="lambda.min")[-1]-b0)^2
err0=mean((Y.test-pred0)^2)^0.5
TPR0=EvaluationMeasures.TPR(Real=trueb, Predicted = as.numeric(round(coef(cvfit,s="lambda.min")[-1],3)!=0))
TNR0=EvaluationMeasures.TNR(Real=trueb, Predicted = as.numeric(round(coef(cvfit,s="lambda.min")[-1],3)!=0))

#obtain the initial MM-Ridge estimator
initial=MMRidge(X,Y,0.05)

D=getD1d(p)
beta0=initial$beta
sig=initial$scale

###compute the adaptive weights
weight.beta=1/abs(beta0)^0.5
weight.theta=1/abs(diff(beta0))^0.5
#limit the weight to avoid extreme values
weight.beta=ifelse(weight.beta>5,5,weight.beta)
weight.theta=ifelse(weight.theta>5, 5, weight.theta)
D.weight=diag(as.vector(weight.theta))%*%D


###choose optimal lambda
test1.cv=raflasso.cv(X, Y, K=5, trim=0.8, type="square", D=D, weight=1, lambda1.min=0.1, lambda2.min=0.1, lambda1.max=0.5, 
                     lambda2.max=0.5, nlambda1=3, nlambda2=3, beta0=beta0)
test2.cv=raflasso.cv(X, Y, K=5, trim=0.8, type="square", D=D.weight, weight=weight.beta, lambda1.min=0.1, lambda2.min=0.1, lambda1.max=0.5, 
                     lambda2.max=0.5, nlambda1=3, nlambda2=3, beta0=beta0)  

test3.cv=raflasso.cv(X, Y, K=5, trim=0.8, type="huber", M=1.345, sig=sig, D=D, weight=1, lambda1.min=0.1, lambda2.min=0.01, lambda1.max=0.3, 
                     lambda2.max=0.03, nlambda1=3, nlambda2=3, beta0=beta0)  
test4.cv=raflasso.cv(X, Y, K=5, trim=0.8, type="huber", M=1.345, sig=sig, D=D.weight, weight=weight.beta, lambda1.min=0.1, lambda2.min=0.01, lambda1.max=0.3, 
                     lambda2.max=0.03, nlambda1=3, nlambda2=3, beta0=beta0)  

test5.cv=raflasso.cv(X, Y, K=5, trim=0.8, type="tukey", M=4.685, sig=sig, D=D, weight=1, lambda1.min=0.1, lambda2.min=0.01, lambda1.max=0.5, 
                     lambda2.max=0.03, nlambda1=3, nlambda2=3, beta0=beta0)  
test6.cv=raflasso.cv(X, Y, K=5, trim=0.8, type="tukey", M=4.685, sig=sig, D=D.weight, weight=weight.beta, lambda1.min=0.1, lambda2.min=0.01, lambda1.max=0.3, 
                     lambda2.max=0.03, nlambda1=3, nlambda2=3, beta0=beta0) 

###obtain results

test1=raflasso(X, Y, type="square", D=D, lambda1=test1.cv$lambda1, lambda2=test1.cv$lambda2, weight=1, beta0=beta0)
err1=mean((Y.test-X.test%*%test1$beta)^2)^0.5
TPR1=EvaluationMeasures.TPR(Real=trueb, Predicted = as.numeric(round(test1$beta,3)!=0))
TNR1=EvaluationMeasures.TNR(Real=trueb, Predicted = as.numeric(round(test1$beta,3)!=0))
rand1=adj.rand.index(group[1:s], cutree(hclust(dist(test1$beta[1:s])),h=h))
norm1=sum((test1$beta-b0)^2)

test2=raflasso(X, Y, type="square", D=D.weight, lambda1=test2.cv$lambda1, lambda2=test2.cv$lambda2, weight=weight.beta, beta0=beta0)
err2=mean((Y.test-X.test%*%test2$beta)^2)^0.5
TPR2=EvaluationMeasures.TPR(Real=trueb, Predicted = as.numeric(round(test2$beta,3)!=0))
TNR2=EvaluationMeasures.TNR(Real=trueb, Predicted = as.numeric(round(test2$beta,3)!=0))
rand2=adj.rand.index(group[1:s], cutree(hclust(dist(test2$beta[1:s])),h=h))
norm2=sum((test2$beta-b0)^2)

test3=raflasso(X, Y, type="huber", M=1.345, sig=sig, D=D, lambda1=test3.cv$lambda1, lambda2=test3.cv$lambda2, weight=1, beta0=beta0)
err3=mean((Y.test-X.test%*%test3$beta)^2)^0.5
TPR3=EvaluationMeasures.TPR(Real=trueb, Predicted = as.numeric(round(test3$beta,3)!=0))
TNR3=EvaluationMeasures.TNR(Real=trueb, Predicted = as.numeric(round(test3$beta,3)!=0))
rand3=adj.rand.index(group[1:s], cutree(hclust(dist(test3$beta[1:s])),h=h))
norm3=sum((test3$beta-b0)^2)

test4=raflasso(X, Y, type="huber", M=1.345, sig=sig, D=D.weight, lambda1=test4.cv$lambda1, lambda2=test4.cv$lambda2, weight=weight.beta, beta0=beta0)
err4=mean((Y.test-X.test%*%test4$beta)^2)^0.5
TPR4=EvaluationMeasures.TPR(Real=trueb, Predicted = as.numeric(round(test4$beta,3)!=0))
TNR4=EvaluationMeasures.TNR(Real=trueb, Predicted = as.numeric(round(test4$beta,3)!=0))
rand4=adj.rand.index(group[1:s], cutree(hclust(dist(test4$beta[1:s])),h=h))
norm4=sum((test4$beta-b0)^2)

test5=raflasso(X, Y, type="tukey", M=4.685, sig=sig, D=D, lambda1=test5.cv$lambda1, lambda2=test5.cv$lambda2, weight=1, beta0=beta0)
err5=mean((Y.test-X.test%*%test5$beta)^2)^0.5
TPR5=EvaluationMeasures.TPR(Real=trueb, Predicted = as.numeric(round(test5$beta,3)!=0))
TNR5=EvaluationMeasures.TNR(Real=trueb, Predicted = as.numeric(round(test5$beta,3)!=0))
rand5=adj.rand.index(group[1:s], cutree(hclust(dist(test5$beta[1:s])),h=h))
norm5=sum((test5$beta-b0)^2)

test6=raflasso(X, Y, type="tukey", M=4.685, sig=sig, D=D.weight, lambda1=test6.cv$lambda1, lambda2=test6.cv$lambda2, weight=weight.beta, beta0=beta0)
err6=mean((Y.test-X.test%*%test6$beta)^2)^0.5
TPR6=EvaluationMeasures.TPR(Real=trueb, Predicted = as.numeric(round(test6$beta,3)!=0))
TNR6=EvaluationMeasures.TNR(Real=trueb, Predicted = as.numeric(round(test6$beta,3)!=0))
rand6=adj.rand.index(group[1:s], cutree(hclust(dist(test6$beta[1:s])),h=h))
norm6=sum((test6$beta-b0)^2)

normball[k,]=c(norm0^0.5,norm1^0.5,norm2^0.5,norm3^0.5,norm4^0.5,norm5^0.5,norm6^0.5)
MSPEall[k,]=c(err0,err1,err2,err3,err4,err5,err6)
TPRall[k,]=c(TPR0,TPR1,TPR2,TPR3,TPR4,TPR5,TPR6)
TNRall[k,]=c(TNR0,TNR1,TNR2,TNR3,TNR4,TNR5,TNR6)
randall[k,]=c(NA,rand1,rand2,rand3,rand4,rand5,rand6)
}


normbave=apply(normball,2,mean)
MSPEave=apply(MSPEall,2,mean)
TPRave=apply(TPRall,2,mean)
TNRave=apply(TNRall,2,mean)
randave=apply(randall,2,mean)

rbind(normbave,TPRave,TNRave,MSPEave,randave)

#write.table(cbind(normball, TPRall, TNRall, MSPEall, randall),file="p50n100xy.txt",append=FALSE,col.names=FALSE,row.names=FALSE)