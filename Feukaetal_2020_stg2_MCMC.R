### 2nd Stage of MCMC Algorithm for Hierarchical State-Space Mvmt Model
### for Multiple Snakes in Parallel
### Update all parameters individually

### Author: Abbey Feuka
### Date: 01JUN20
### Notes: 
stsp.var.2.indv <- function(Delta,s20,g.mat,t.mat,b.mat,
                            s21.mat,X,n.mcmc){ 
  n=dim(g.mat)[1] #number of individuals
  n.mcmc1=dim(g.mat)[2] #number of mcmc from stage 1
  
  ####
  ####  Libraries
  ####
  library(msm)
  library(MCMCpack)
  library(CircStats)
  library(parallel)
  library(mvtnorm)
  library(mvnfast)
  library(boot)
  
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
  logit=function(x){
    y=log(x/(1-x))
    return(y)
  }
  dIG <- function(x,q,r){
    x^(-q-1) * exp(-1/r/x) / (r^q) / gamma(q)
  }
  
  llh<-function(Delta,gamma,theta,s21,p){
    x=0
    M=make.M(gamma,theta)
    for(t in 2:nrow(Delta)){
      tmp.diff=Delta[t,]-M%*%Delta[t-1,]
      x=x+log(p[t]*sum(dnorm(tmp.diff,0,sqrt(s21)))+
                    (1-p[t])*sum(dnorm(Delta[t,],0,sqrt(s20))))
    }
    return(x)
  }
  ####
  ####  Set Up Variables
  ####
  pp=dim(b.mat)[2] #number of coefficients
  beta.save=array(NA,dim=c(n.mcmc,pp,n))
  p.save=lapply(1:n,function(i){matrix(NA,n.mcmc,nrow(na.omit(X[,,i])))})
  theta.save=gamma.save=s21.save=matrix(NA,n,n.mcmc)
  mu.b.save=matrix(NA,n.mcmc,pp)
  Sigma.b.save=array(NA,dim=c(pp,pp,n.mcmc))
  mu.g.save=s2.g.save=mu.t.save=s2.t.save=
    mu.logs1.save=s2.logs1.save=rep(NA,n.mcmc)

  ##Stage 1 Priors 
  mu.b.1st=rep(0,pp)
  Sigma.b.1st=1*diag(pp)
  Sigma.b.1st.inv=solve(Sigma.b.1st)
  q1.1st=1
  r1.1st=0.1
  mu.t.1st=0
  rho.t.1st=0.1
  alpha.g1.1st=1
  alpha.g2.1st=1
  
  ####
  ####  Starting Values
  ####
  ###Indiv params (param_ij)
  beta=beta.tmp=t(sapply(1:n,function(i){rep(0,pp)}))
  p=p.tmp=lapply(1:n,function(i){c(NA,pnorm(na.omit(X[,,i])%*%beta[i,]))})
  theta=rep(0,n)
  gamma=rep(0.8,n)
  s21=rep(100,n)

  ###Treatment params (param_j)
  mu.b=rep(0,pp) #int, trees, building, road, grass
  Sigma.b=1*diag(pp)
  Sigma.b.inv=solve(Sigma.b)
  mu.g=logit(0.8)
  s2.g=1.25
  mu.t=pi
  s2.t=(pi^2)/4
  mu.logs1=log(20) #log(SD)
  s2.logs1=1
  
  ###
  ### Hyperpriors (param_pop)
  ###
  #means
  mu.b.pop=rep(0,pp) 
  Sigma.b.p=1*diag(pp)
  Sigma.b.p.inv=solve(Sigma.b.p)
  mu.g.pop=logit(0.7)
  s2.g.pop=1.5
  mu.t0=pi
  s2.t0=(pi^2)/4
  mu.mus1.pop=log(20)
  s2.mus1.pop=3
  
  #variances
  S=matrix(0,pp,pp)
  diag(S)=1
  nu=3 
  q.g=0.01
  r.g=10
  alpha.rho1=1
  alpha.rho2=1
  q.s21.pop=0.0001
  r.s21.pop=1000
  q.t=2.5
  r.t=0.5
  
  ####
  ####  Begin MCMC Loop
  ####
  for(k in 1:n.mcmc){
    if(k%%10000==0) cat(k," ")
    
    ####
    ####  Individual-level parameters
    ####  MH Updates Propose from 1st Stage Samples
    #beta
    for(i in 1:n){
    #loop over coefficients
      for(j in 1:pp){
        nj=seq(1,pp)[-j]
        beta.star=beta.tmp[i,]
        samp=sample(1:(n.mcmc1-1),1,replace=T)
        beta.star[j]=b.mat[samp,j,i]
        p.star=c(NA,pnorm(na.omit(X[,,i])%*%beta.star))
        Sig.uo=t(as.matrix(Sigma.b[-nj,nj]))
        Sig.oo.inv=solve(Sigma.b[nj,nj])
        Sig.uo.1st=t(as.matrix(Sigma.b.1st[-nj,nj]))
        Sig.oo.inv.1st=solve(Sigma.b.1st[nj,nj])
        
        mh1=llh(Delta[[i]],gamma[i],theta[i],s21[i],p.star)+
          dnorm(beta.star[j],
                mu.b[j]+(Sig.uo)%*%Sig.oo.inv%*%(beta.tmp[i,nj]-mu.b[nj]),
                sqrt(Sigma.b[j,j]-(Sig.uo)%*%Sig.oo.inv%*%t(Sig.uo)),log=T)+
          llh(Delta[[i]],gamma[i],theta[i],s21[i],p[[i]])+
          dnorm(beta[i,j],
                mu.b.1st[j]+(Sig.uo.1st)%*%Sig.oo.inv.1st%*%(beta.tmp[i,nj]-mu.b.1st[nj]), #beta or beta.tmp?
                sqrt(Sigma.b.1st[j,j]-(Sig.uo.1st)%*%Sig.oo.inv.1st%*%t(Sig.uo.1st)),log=T)
        mh2=llh(Delta[[i]],gamma[i],theta[i],s21[i],p.tmp[[i]])+
          dnorm(beta[i,j],
                mu.b[j]+(Sig.uo)%*%Sig.oo.inv%*%(beta.tmp[i,nj]-mu.b[nj]),
                sqrt(Sigma.b[j,j]-(Sig.uo)%*%Sig.oo.inv%*%t(Sig.uo)),log=T)+
          llh(Delta[[i]],gamma[i],theta[i],s21[i],p.star)+
          dnorm(beta.star[j],
                mu.b.1st[j]+(Sig.uo.1st)%*%Sig.oo.inv.1st%*%(beta.star[nj]-mu.b.1st[nj]),
                sqrt(Sigma.b.1st[j,j]-(Sig.uo.1st)%*%Sig.oo.inv.1st%*%t(Sig.uo.1st)),log=T)
        
        mh=exp(mh1-mh2)
        if(mh>runif(1) & !is.na(mh)){
          beta.tmp[i,j]=beta.star[j]
          p.tmp[[i]]=p.star
        }
      }
      beta[i,]=beta.tmp[i,]
      p[[i]]=p.tmp[[i]]

    #gamma
    gamma.star=g.mat[i,sample(n.mcmc1,1,replace=T)]
    mh1=llh(Delta[[i]],gamma.star,theta[i],s21[i],p[[i]])+
      log(dnorm(logit(gamma.star),mu.g,sqrt(s2.g))* #[logit(gamma)*_ij|mu.g_j,s2.g_j]
          dbeta(gamma[i],alpha.g1.1st,alpha.g2.1st)*
          exp(logit(gamma[i]))/(1+exp(logit(gamma[i])))^2) #|d/dlogit(gamma) g(logit(gamma))|
    mh2=llh(Delta[[i]],gamma[i],theta[i],s21[i],p[[i]])+
      log(dnorm(logit(gamma[i]),mu.g,sqrt(s2.g))* #[logit(gamma)*_ij|mu.g_j,s2.g_j]
          dbeta(gamma.star,alpha.g1.1st,alpha.g2.1st)*
          exp(logit(gamma.star))/(1+exp(logit(gamma.star)))^2) #|d/dlogit(gamma) g(logit(gamma))|
    mh=exp(mh1-mh2)
    if(mh>runif(1) & !is.na(mh)){
      gamma[i]=gamma.star
    }
    
    #theta
    theta.star=t.mat[i,sample(n.mcmc1,1,replace=T)]
    mh1=llh(Delta[[i]],gamma[i],theta.star,s21[i],p[[i]])+
      log(dnorm(theta.star,mu.t,sqrt(s2.t))*
          dwrpcauchy(theta[i],mu.t.1st,rho.t.1st)) #1st stage prior
    mh2=llh(Delta[[i]],gamma[i],theta[i],s21[i],p[[i]])+
      log(dnorm(theta[i],mu.t,sqrt(s2.t))*
          dwrpcauchy(theta.star,mu.t.1st,rho.t.1st)) #1st stage prior
    mh=exp(mh1-mh2)
    if(mh>runif(1) & !is.na(mh)){
      theta[i]=theta.star
    }
    
    #log(sigma1)
    s21.star=s21.mat[i,sample(n.mcmc1,1,replace=T)]
    mh1=llh(Delta[[i]],gamma[i],theta[i],s21.star,p[[i]])+
      log(dnorm(log(sqrt(s21.star)),mu.logs1,sqrt(s2.logs1))*
          dIG(s21[i],q=q1.1st,r=r1.1st)* #1st stage prior
          s21[i]) #|d/dlog(s21) g(log(s21)|
    mh2=llh(Delta[[i]],gamma[i],theta[i],s21[i],p[[i]])+
      log(dnorm(log(sqrt(s21[i])),mu.logs1,sqrt(s2.logs1))*
            dIG(s21.star,q=q1.1st,r=r1.1st)* #1st stage prior
            s21.star) #|d/dlog(s21) g(log(s21)|
    mh=exp(mh1-mh2)
    if(mh>runif(1) & !is.na(mh)){
      s21[i]=s21.star
    }
    } #end of indiv loop

    ####
    #### Sample mu.b_j
    #### Gibbs
    Sigma.tmp=solve(n*Sigma.b.inv+Sigma.b.p.inv) #A^-1
    b.sum=colSums(beta)
    mu.tmp=Sigma.tmp%*%t(b.sum%*%Sigma.b.inv+mu.b.pop%*%Sigma.b.p.inv) #A^-1 * b
    mu.b=as.vector(rmvn(1,mu.tmp,Sigma.tmp))
    
    ####
    #### Sample Sigma.b_j
    #### Gibbs
    S.b=matrix(0,pp,pp)
    for(i in 1:n){
      S.b= S.b+(beta[i,]-mu.b)%*%t(beta[i,]-mu.b)
    }
    S.tmp=solve(S.b+S*nu)
    nu.tmp=nu+n
    Sigma.b.inv=rwish(nu.tmp,S.tmp)
    Sigma.b=solve(Sigma.b.inv)
    
    ####
    #### Sample mu.g_j
    #### Gibbs
    tmp.sum=sum(logit(gamma))/s2.g
    tmp.var=1/(n/s2.g+1/s2.g.pop) #A^-1
    tmp.mn=tmp.var*(tmp.sum+ mu.g.pop/s2.g.pop)#A^-1*b
    mu.g=rnorm(1,tmp.mn,sqrt(tmp.var)) #on logit scale
    
    ###
    ### Sample s2.g_j
    ### Gibbs
    tmp.q=n/2+q.g
    tmp.r=sum((logit(gamma)-mu.g)^2)/2 + 1/r.g
    s2.g=1/rgamma(1,tmp.q,,1/tmp.r)
    
    ####
    #### Sample mu.t_j
    #### Gibbs
    tmp.var=1/(n/s2.t+1/s2.t0)
    tmp.mn=tmp.var*(sum(theta)/s2.t+mu.t0/s2.t0)
    mu.t=rnorm(1,tmp.mn,sqrt(tmp.var))
    
    ####
    #### Sample s2.t_j
    #### Gibbs
    tmp.q=n/2+q.t
    tmp.r=1/(sum((theta-mu.t)^2)/2+1/r.t)
    s2.t=1/rgamma(1,tmp.q,,tmp.r)
    
    ###
    ### Sample mu.log(sigma1_j)
    ### Gibbs
    tmp.sum=sum(log(sqrt(s21)))/s2.logs1
    tmp.var=1/(n/s2.logs1+1/s2.mus1.pop) #A^-1
    tmp.mn=tmp.var*(tmp.sum+mu.mus1.pop/s2.mus1.pop)#A^-1*b
    mu.logs1=rnorm(1,tmp.mn,sqrt(tmp.var))
    
    ###
    ### Sample s2.log(sigma1_j)
    ### Gibbs
    tmp.q=n/2+q.s21.pop
    tmp.r=sum((log(sqrt(s21))-mu.logs1)^2)/2 + 1/r.s21.pop
    s2.logs1=1/rgamma(1,tmp.q,,1/tmp.r)
    
    ####
    ####  Save Samples
    ####
      theta.save[,k]=theta
      gamma.save[,k]=gamma
      s21.save[,k]=s21
      mu.t.save[k]=mu.t
      s2.t.save[k]=s2.t
      mu.g.save[k]=mu.g
      s2.g.save[k]=s2.g 
      mu.logs1.save[k]=mu.logs1
      s2.logs1.save[k]=s2.logs1 
      for(i in 1:n){
        beta.save[k,,i]=beta[i,]
        p.save[[i]][k,]=p[[i]][-1]
      }
      mu.b.save[k,]=mu.b
      Sigma.b.save[,,k]=Sigma.b

  }
  cat("\n")
  
  ####
  #### Write Output
  ####
  
  list(s21.save=s21.save,theta.save=theta.save,
       gamma.save=gamma.save,beta.save=beta.save,p.save=p.save,
       mu.b.save=mu.b.save,Sigma.b.save=Sigma.b.save,
       mu.t.save=mu.t.save,s2.t.save=s2.t.save,
       mu.g.save=mu.g.save,s2.g.save=s2.g.save,
       mu.logs1.save=mu.logs1.save,s2.logs1.save=s2.logs1.save,
       n.mcmc=n.mcmc)
}