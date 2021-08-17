

#setwd("C:/Users/carlos.maza/Desktop/App")
library(ccgarch)
library(tseries) 
library(data.table)
source("/Users/ctruciosm/Desktop/GDFM-CHF/aux_functions/dcc_estimation2.R")
source("/Users/ctruciosm/Desktop/GDFM-CHF/aux_functions/loglik_dcc2.R")
source("/Users/ctruciosm/Desktop/GDFM-CHF/aux_functions/R_dcc.R")
source("/Users/ctruciosm/Desktop/GDFM-CHF/aux_functions/grad_dcc2.R")
source("/Users/ctruciosm/Desktop/GDFM-CHF/aux_functions/PVC_Qiwei.R")
source("/Users/ctruciosm/Desktop/GDFM-CHF/aux_functions/NumComp.R")
# Estimator ratio 1
# Kaisse-Gutman 32
# Bai (8 when the maximun is 8), 656 when maxiumum 656, se maximun 640 then 28 VC.
retornos = fread("/Users/ctruciosm/Desktop/GDFM-CHF/Data/retornos_APP3.txt")[,-1]
#Y = scale(as.matrix(retornos),scale=FALSE)
#m1y = comVolQ(Y, m = 5)
#NumComp(m1y$values/sum(m1y$values))

#M1 = m1y$vectors
#Fa = Y%*%M1
#N = dim(retornos)[2]
#Tt = dim(retornos)[1]
#Bai = c()
#for (i in 1:N){
#   Bai[i] = log(mean((Y-t(M1[,1:i]%*%t(Fa[,1:i])))^2)) + i*log(min(N,Tt))/(N*Tt/(N+Tt))
# }
#which(Bai[1:640] == min(Bai[1:640]))


n=dim(retornos)[1]; p=dim(retornos)[2]
#Setting parameters for the RW scheme.
W = 750; L = n;  WR = L-750
Hvector = matrix(NA,nrow=p*p,ncol=WR)
p = 8
partes = 1:50
for (l in partes){
  print(l)
  H = matrix(0,p,p)
  res = matrix(0,ncol=p,nrow=750)
  h = c()
  Y = scale(as.matrix(retornos[l:(l+750-1),]),scale=FALSE)
  m1y = comVolQ(Y, m = 5)
  M1 = m1y$vectors[,1:p]
  M2 = m1y$vectors[,-c(1:p)]
  X = Y%*%M1 #Volatility Components
  X_c = X #-mean(X)

  for (i in 1:p){
  GARCHaux1 = garch(X[,i], order = c(1,1),trace=FALSE)
  if(sum(GARCHaux1$coef[2:3])<1){
    H[i,i] = sqrt(GARCHaux1$coef[1] + GARCHaux1$coef[2]*X[750,i]^2+GARCHaux1$coef[3]*GARCHaux1$fitted.values[750,1]^2)
    res[,i] = c(as.vector(GARCHaux1$residuals[-1]),0)
  }else{
    param = GARCHestim(X[,i])
    h[1] = var(X[,i])
    for (ii in 2:751){
      h[ii] = param[1]+param[2]*X[ii-1,i]^2+param[3]*h[ii-1]
    }
    H[i,i] = sqrt(h[751])
    res[,i] = c(X[2:750,i]/sqrt(h[2:750]),0)
  }
  }
  
  parDCC1 = dcc.estimation2(res[1:749,], para = c(0.05,0.9), gradient = 1)$par
  if(sum(parDCC1)>1){
    parDCC1 = dcc.estimation2(stdres1[1:749,], para = c(0.1,0.1), gradient = 1)$par
  }
  uncR <- cov(res[1:749,])
  out <- .Call("dcc_est", res, uncR, parDCC1[1], parDCC1[2])
  R1 = out[[1]]
  H[1:p,1:p] = H[1:p,1:p]%*%matrix(R1[750,],p,p)%*%H[1:p,1:p]
  
  
  Hhat = M1%*%H%*%t(M1) + m1y$COV%*%M2%*%t(M2) + M2%*%t(M2)%*%m1y$COV%*%M1%*%t(M1)  
  Hvector[,l] = as.vector(Hhat)

}

Sigma = Hvector[,partes]
write.table(data.frame(Sigma),"H_GPVC_P1.txt",sep=",",dec=".",col.names=FALSE,row.names = FALSE)


