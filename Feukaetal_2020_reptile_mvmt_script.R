library(sf)
library(rgdal)
library(raster)
library(dplyr)
library(parallel)
library(circular)
library(boot)

setwd("~/OneDrive - Colostate/Manuscripts/EA 2020/Final Code")

###
### Movement Data
###
#movement data
bts = read.csv("Feukaetal_2020_reptile_mvmt.csv")
n <- length(unique(bts$idx))

#convert to lists of dates and deltas
bts.delta = list()
bts.date = list()
for(i in 1:n){
  bts.delta[[i]]=cbind.data.frame(bts$Delta_easting[bts$idx==i],bts$Delta_northing[bts$idx==i])
  bts.date[[i]]=as.Date(bts$date[bts$idx==i],format="%m/%d/%y")
}

#trim to just covariates
X=list()
for(i in 1:n){
  X[[i]]=bts[bts$idx==i,8:11] #columns with covariates
}

#line up with each step's starting position (t-1)
for(i in 1:n){
  X[[i]]=X[[i]][1:(nrow(X[[i]])-1),1:4]
}

#lists to arrays
x.test=numeric(n)
for(i in 1:n){
  x.test[i]=nrow(X[[i]])
}

X.a=array(NA,dim=c(max(x.test),dim(X[[1]])[2],n)) #discrete points
for(i in 1:n){
  X.a[1:nrow(X[[i]]),,i]=as.matrix(X[[i]])
}

#############
# 1st Stage
#############

#individuals that did not move enough to be classified as in movement stage
#manually set these individuals to always be in residence state
no.mov=as.factor(c(4016,4017,4018,4031,4062,4064,4074,5007,5018,5021,5023))
no.mov.id=as.factor(sort(as.character(no.mov)))

#fit first stage models in parallel
n.mcmc1=100000
source('Feukaetal_2020_stg1_MCMC.R')
stg1=mclapply(1:n,function(i){
  stsp.var.1.z0(Delta=bts.delta[[i]],s20=25,X=X[[i]],
                n.mcmc=n.mcmc1,no.mov=(i %in% no.mov.id))},
  mc.cores = 8)

#save first stage samples
{gamma_ij=theta_ij=s21_ij=matrix(0,length(bts.delta),n.mcmc1)
  beta_ij=array(0,dim=c(n.mcmc1,dim(X[[1]])[2],n))
  Delta=list()
  for(i in 1:n){
    gamma_ij[i,]=stg1[[i]]$gamma.save
    theta_ij[i,]=stg1[[i]]$theta.save
    s21_ij[i,]=stg1[[i]]$s21.save
    beta_ij[,1:4,i]=stg1[[i]]$beta.save[,1:4]
    Delta[[i]]=stg1[[i]]$Delta 
  }
}

#############
# 2nd Stage
#############
t.id=aggregate(bts$idx, by=list(bts$tidx,bts$idx), FUN=max) 
trmt=t.id$Group.1 #vector of treatment indicies CU=1, FU=2, UU=3

#separate movements into list
Delta1=Delta2=Delta3=list()
for(i in 1:n){
  if(trmt[i]==1){
    Delta1[[i]]=bts.delta[[i]]
  }
}
for(i in 1:n){
  if(trmt[i]==2){
    Delta2[[i]]=bts.delta[[i]]
  }
}
for(i in 1:n){
  if(trmt[i]==3){
    Delta3[[i]]=bts.delta[[i]]
  }
}

names(Delta1) <- seq_along(Delta1)
Delta1[sapply(Delta1,is.null)] <- NULL
length(Delta1)
names(Delta2) <- seq_along(Delta2)
Delta2[sapply(Delta2,is.null)] <- NULL
length(Delta2)
names(Delta3) <- seq_along(Delta3)
Delta3[sapply(Delta3,is.null)] <- NULL
length(Delta3)

Delta.t=list(Delta1,Delta2,Delta3)

#fit second stage models
n.mcmc2=50000
source('Feukaetal_2020_stg2_MCMC')
test2=mclapply(1:length(unique(trmt)),function(i){
  set.seed(1)
  stsp.var.2.indv(Delta=Delta.t[[i]],s20=25,g.mat=gamma_ij[trmt==i,],
                  t.mat=theta_ij[trmt==i,],b.mat=beta_ij[,,trmt==i],
                  s21.mat=s21_ij[trmt==i,],X=X.a[,,trmt==i],
                  n.impt=5000,n.mcmc=n.mcmc2)
},mc.cores=3)

#save second stage samples
gamma2=theta2=s21.2=matrix(0,length(bts.delta),n.mcmc2)
mu.g2=mu.t2=mu.logs1.2=s2.g2=s2.t2=s2.logs1.2=matrix(0,length(unique(trmt)),n.mcmc2)
mu.beta2=array(0,dim=c(n.mcmc2,dim(beta_ij)[2],length(unique(trmt))))
for(j in 1:length(unique(trmt))){
  gamma2[trmt==j,]=test2[[j]]$gamma.save
  theta2[trmt==j,]=test2[[j]]$theta.save
  s21.2[trmt==j,]=test2[[j]]$s21.save
  mu.g2[j,]=test2[[j]]$mu.g.save
  mu.t2[j,]=test2[[j]]$mu.t.save
  mu.logs1.2[j,]=test2[[j]]$mu.logs1.save
  mu.beta2[,1:dim(beta_ij)[2],j]=test2[[j]]$mu.b.save[,1:dim(beta_ij)[2]]
  s2.g2[j,]=test2[[j]]$s2.g.save
  s2.t2[j,]=test2[[j]]$s2.t.save
  s2.logs1.2[j,]=test2[[j]]$s2.logs1.save
}
pp=dim(beta_ij)[2]
beta2=array(NA,dim=c(n.mcmc2,pp,length(bts.delta)))
for(j in 1:length(unique(trmt))){
  beta2[,1:pp,trmt==j]=test2[[j]]$beta.save[,1:pp,]
}

##############
# Remove Burn in
##############
#1st stage
n.burn=n.mcmc1/10
gamma1=theta1=s21.1=matrix(NA,length(bts.delta),n.mcmc1-n.burn)
beta1=array(NA,dim=c(n.mcmc1-n.burn,dim(X[[1]])[2],length(bts.delta)))
for(i in 1:length(bts.delta)){
  gamma1[i,]=gamma_ij[i,-(1:n.burn)]
  theta1[i,]=theta_ij[i,-(1:n.burn)]
  s21.1[i,]=s21_ij[i,-(1:n.burn)]
  beta1[,1:dim(X[[1]])[2],i]=beta_ij[-(1:n.burn),1:dim(X[[1]])[2],i]
}

#second stage
n.burn2=n.mcmc2/10
gamma2=theta2=s21.2=matrix(0,length(bts.delta),n.mcmc2-n.burn2)
mu.g2=mu.t2=mu.logs1.2=s2.g2=s2.t2=s2.logs1.2=matrix(0,length(unique(trmt)),n.mcmc2-n.burn2)
mu.beta2=array(0,dim=c(n.mcmc2-n.burn2,dim(beta_ij)[2],length(unique(trmt))))
for(j in 1:length(unique(trmt))){
  gamma2[trmt==j,]=test2[[j]]$gamma.save[,-(1:n.burn2)]
  theta2[trmt==j,]=test2[[j]]$theta.save[,-(1:n.burn2)]
  s21.2[trmt==j,]=test2[[j]]$s21.save[,-(1:n.burn2)]
  mu.g2[j,]=test2[[j]]$mu.g.save[-(1:n.burn2)]
  mu.t2[j,]=test2[[j]]$mu.t.save[-(1:n.burn2)]
  mu.logs1.2[j,]=test2[[j]]$mu.logs1.save[-(1:n.burn2)]
  mu.beta2[,1:dim(beta_ij)[2],j]=test2[[j]]$mu.b.save[-(1:n.burn2),1:dim(beta_ij)[2]]
  s2.g2[j,]=test2[[j]]$s2.g.save[-(1:n.burn2)]
  s2.t2[j,]=test2[[j]]$s2.t.save[-(1:n.burn2)]
  s2.logs1.2[j,]=test2[[j]]$s2.logs1.save[-(1:n.burn2)]
}

#transforming vars using MC sampling
t.mu.g2=t.mu.s2=t.s2.g2=t.s2.s2=t.mu.s=t.s2.s=matrix(0,length(unique(trmt)),length(mu.g2[1,]))
for(j in 1:length(unique(trmt))){
  for(k in 1:length(mu.g2[1,])){
    t.mu.g2[j,k]=mean(inv.logit(rnorm(10000,mu.g2[j,k],sqrt(s2.g2[j,k]))))
    t.s2.g2[j,k]=var(inv.logit(rnorm(10000,mu.g2[j,k],sqrt(s2.g2[j,k]))))
    t.mu.s2[j,k]=mean(exp(2*rnorm(10000,mu.logs1.2[j,k],sqrt(s2.logs1.2[j,k]))))
    t.s2.s2[j,k]=var(exp(2*rnorm(10000,mu.logs1.2[j,k],sqrt(s2.logs1.2[j,k]))))
    t.mu.s[j,k]=mean(exp(rnorm(10000,mu.logs1.2[j,k],sqrt(s2.logs1.2[j,k]))))
    t.s2.s[j,k]=var(exp(rnorm(10000,mu.logs1.2[j,k],sqrt(s2.logs1.2[j,k]))))
  }
}

#law of total variance
var.s=numeric(3)
for(j in 1:length(unique(trmt))){
  var.s[j-1]=mean(t.s2.s[j,])+var(t.mu.s[j,])
}

t.mu.p=t.s2.p=array(NA,dim=dim(mu.beta2))
for(j in 1:length(unique(trmt))){
  Sig=test2[[j]]$Sigma.b.save[,,-(1:n.burn2)]
  for(k in 1:dim(mu.beta2)[1]){
    t.mu.p[k,1:4,j] <- colMeans(pnorm(rmvn(10000,mu.beta2[k,1:4,j],Sig[,,k])))
    t.s2.p[k,1:4,j] <- apply(pnorm(rmvn(10000,mu.beta2[k,1:4,j],Sig[,,k])),2,var)
  }
}

pp=dim(beta_ij)[2]
beta2=array(NA,dim=c(n.mcmc2-n.burn2,pp,length(bts.delta)))
for(j in 1:length(unique(trmt))){
  beta2[,1:pp,trmt==j]=test2[[j]]$beta.save[-(1:n.burn2),1:pp,]
}
