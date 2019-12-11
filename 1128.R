rm(list=ls())
library(np)
data('cps71')
head(cps71)
#1
plot(cps71$age, cps71$lowage, 
     xlab = 'Age', ylab = 'lowage')

#2
regressogram = function(x,y,left,right,k,plotit,
                        xlab='',ylab=''){
  n = length(x)
  B = seq(left, right, length = k+1)
  WhichBin = findInterval(x,B)
  N = tabulate(WhichBin)
  m.hat = rep(0,k)
  for(j in 1:k){
    if(N[j]>0)m.hat[j]=mean(y[WhichBin == j])
  }
  if(plotit==TRUE){
    a = min(c(y,m.hat))
    b = max(c(y,m.hat))
    plot(B,c(m.hat, m.hat[k]),lwd=3, type ='s',
         xlab=xlab, ylab=ylab, ylim=c(a,b),
         col='blue')
    points(x,y)
  }
  return(list(bins=B,m.hat=m.hat))
}

age <- cps71$age
lowage <- cps71$logwage
out = regressogram(age,lowage,left=21,right=65,k=20,plotit=TRUE)

#3
kernel = function(x,y,grid,h){
  n = length(x)
  k = length(grid)
  m.hat = rep(0,k)
  for(i in 1:k){
    w = dnorm(grid[i],x,h)
    m.hat[i] = sum(y*w)/sum(w)
  }
  return(m.hat)
}
kernel.fitted = function(x,y,h){
  n = length(x)
  m.hat = rep(0,n)
  S = rep(0,n)
  for(i in 1:n){
    w = dnorm(x[i],x,h)
    w = w/sum(w)
    m.hat[i] = sum(y*w)
    S[i] = w[i]
  }
  return(list(fitted=m.hat,S=S))
}

CV = function(x,y,H){
  n = length(x)
  k = length(H)
  cv = rep(0,k)
  nu = rep(0,k)
  gcv = rep(0,k)
  for(i in 1:k){
    tmp = kernel.fitted(x,y,H[i])
    cv[i] = mean(((y - tmp$fitted)/(1-tmp$S))^2)
    nu[i] = sum(tmp$S)
    gcv[i] = mean((y - tmp$fitted)^2)/(1-nu[i]/n)^2
  }
  return(list(cv=cv,gcv=gcv,nu=nu))
}
H = seq(.1,5,length=20)
out = CV(age, lowage,H)
j = which.min(out$cv)
h = H[j]
h

#4
h
grid = seq(min(age),max(age),length=200)
m.hat = kernel(age, lowage,grid,h)
plot(age, lowage ,xlab="Age",ylab="lowage)")
lines(grid,m.hat,lwd=3,col="blue")
lines(h,CV$gcv,lty=2,col="red",lwd=3)
lines(CV$nu,CV$gcv,lty=2,col="red",lwd=3)
