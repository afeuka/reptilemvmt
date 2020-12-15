###
### 1st Stage MCMC Algorithm for State-Space Mvmt Model for One Snake
### Using Gibbs Updates on Gamma (unif prior), s20 and s21 (IG priors), importance sampling on theta and beta
### Lunn Method 
### Fit this model in parallel individually for each snake

### Author: Abbey Feuka
### Date: 15JUN20
### Notes: Modification of 07APR20 script
### Filters out low-movement indivs and sets all z==0 (keeps s2_1>s2_0)

stsp.var.1.z0 <- function(Delta,s20,X,n.mcmc,no.mov=F){
  
  ####
  ####  Libraries
  ####
  library(msm)
  library(CircStats)
  library(mvnfast)
  
  make.Theta <- function(theta){
    Theta=matrix(0,2,2) 
    Theta[1,1]=cos(theta)
    Theta[1,2]=-sin(theta)
    Theta[2,1]=sin(theta)
    Theta[2,2]=cos(theta)
    Theta 
  }
  make.M <- function(gamma,theta){
    Theta=make.Theta(theta) 
    M=gamma*Theta 
    M
  }
  
  ####
  ####  Set Up Variables
  ####
  
  X=as.matrix(X)
  pp=dim(X)[2] #number of coefficients
  Delta=as.matrix(Delta)
  Td=dim(Delta)[1] #timesteps
  
  theta.save=rep(0,n.mcmc)
  gamma.save=rep(0,n.mcmc)
  s21.save=rep(0,n.mcmc)
  beta.save=matrix(0,n.mcmc,pp)
  z.save=p.save=matrix(NA,Td,n.mcmc)
  
  ####
  ####  Priors and Starting Values
  ####
  #starting values
  theta=0
  gamma=0.9
  Theta=make.Theta(theta)
  M=make.M(gamma,theta)
  s21=10
  Sig1=s21*diag(2)
  Sig1.inv=solve(Sig1)
  beta=rep(0,pp)
  p=c(1,pnorm(X%*%beta))
  v=rnorm(Td,0,1)
  if(no.mov==FALSE){
    z=rbinom(Td,1,p)
  } else if(no.mov==TRUE){
    z=rep(0,Td)
  }
  
  #hyperpriors for beta
  mu.beta=rep(0,pp)
  Sigma.beta=matrix(0,pp,pp)
  diag(Sigma.beta)=1
  Sig.beta.inv=solve(Sigma.beta)
  q1=1
  r1=0.1
  mu.t=0
  rho.t=0.1
  alpha.g1=1
  alpha.g2=1
  
  #tuning
  t.acc=0
  g.acc=0
  
  ####
  ####  Begin MCMC Loop
  ####
  for(k in 1:n.mcmc){
    if(k%%1000==0) {
      cat(k," ")
    }
    ###
    ### Sample z  
    ### Gibbs 
    if(no.mov==FALSE){
      for(t in 2:Td){
        tmp.diff=Delta[t,]-M%*%Delta[t-1,]
        p.tmp=p[t]*sum(dnorm(tmp.diff,0,sqrt(s21)))/
          (p[t]*sum(dnorm(tmp.diff,0,sqrt(s21)))+(1-p[t])*sum(dnorm(Delta[t,],0,sqrt(s20))))
        z[t]=rbinom(1,1,p.tmp)
      }
    }
    
    ###
    ### Sample v_t
    ### Gibbs
    v[-1][z[-1]==0]=rtnorm(sum(1-z[-1]),(X%*%beta)[z[-1]==0],1,-Inf,0) 
    v[-1][z[-1]==1]=rtnorm(sum(z[-1]),(X%*%beta)[z[-1]==1],1,0,Inf)
    
    ###
    ### Sample beta
    ### Gibbs
    tmp.b=t(v[-1])%*%X+t(mu.beta)%*%Sig.beta.inv
    tmp.a=t(X)%*%X+Sig.beta.inv
    tmp.var=solve(tmp.a)
    tmp.mn=tmp.var%*%t(tmp.b)
    beta=as.vector(rmvn(1,tmp.mn,tmp.var))
    p=c(1,pnorm(X%*%beta))
    
    ####
    ####  Sample s21
    ####  Gibbs
    tmp.diff=Delta[2:Td,]-t(M%*%t(Delta[1:(Td-1),]))
    tmp.sum=sum(tmp.diff[z[-1]==1,]^2)
    s21=1/rgamma(1,sum(z[-1])+q1,,1/(tmp.sum/2+1/r1))
    Sig1=s21*diag(2)
    Sig1.inv=solve(Sig1)
    
    ####
    ####  Sample theta 
    ####  MH Propose from Prior
    theta.star=rwrpcauchy(1,mu.t,rho.t)
    M.star=make.M(gamma,theta.star)
    tmp.theta.sum.1=0
    tmp.theta.sum.2=0
    for(t in 2:Td){
      if(z[t]==1){
        delta.star=M.star%*%Delta[t-1,]
        delta=M%*%Delta[t-1,]
        tmp.theta.sum.1=tmp.theta.sum.1-t(Delta[t,]-delta.star)%*%Sig1.inv%*%(Delta[t,]-delta.star)/2
        tmp.theta.sum.2=tmp.theta.sum.2-t(Delta[t,]-delta)%*%Sig1.inv%*%(Delta[t,]-delta)/2
      }
    }
    mh1=tmp.theta.sum.1
    mh2=tmp.theta.sum.2
    mh=exp(mh1-mh2)
    if(mh > runif(1)){
      theta=theta.star
      Theta=make.Theta(theta.star)
      M=M.star
      t.acc=t.acc+1
    }
    
    ####
    ####  Sample gamma 
    ####  MH Propose from Prior
    gamma.star=rbeta(1,alpha.g1,alpha.g2)
    M.star=make.M(gamma.star,theta)
    tmp.gamma.sum.1=0
    tmp.gamma.sum.2=0
    for(t in 2:Td){
      if(z[t]==1){
        delta.star=M.star%*%Delta[t-1,]
        delta=M%*%Delta[t-1,]
        tmp.gamma.sum.1=tmp.gamma.sum.1-t(Delta[t,]-delta.star)%*%Sig1.inv%*%(Delta[t,]-delta.star)/2
        tmp.gamma.sum.2=tmp.gamma.sum.2-t(Delta[t,]-delta)%*%Sig1.inv%*%(Delta[t,]-delta)/2
      }
    }
    mh1=tmp.gamma.sum.1
    mh2=tmp.gamma.sum.2
    mh=exp(mh1-mh2)
    if(mh > runif(1)){
      gamma=gamma.star
      M=M.star
      g.acc=g.acc+1
    }
    
    ####
    ####  Save Samples
    ####
    theta.save[k]=theta
    gamma.save[k]=gamma
    s21.save[k]=s21
    beta.save[k,]=beta
    z.save[,k]=z
  }
  cat("\n")
  
  #Acceptance
  t.acc=t.acc/n.mcmc
  g.acc=g.acc/n.mcmc
  
  ####
  #### Write Output
  ####
  list(gamma.save=gamma.save,theta.save=theta.save,beta.save=beta.save,
       s21.save=s21.save,z.save=z.save,Delta=Delta,n.mcmc=n.mcmc,
       t.acc=t.acc,g.acc=g.acc)
}
