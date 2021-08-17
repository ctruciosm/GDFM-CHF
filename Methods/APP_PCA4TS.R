#########################################################
#####   Segment Multivriate Volatility Processes  #######
#########################################################
# Loading packages
require(tseries) 
library(PCA4TS)
library(ccgarch) # deprecated
library(ccgarch2)
library(data.table)
source("/Users/ctruciosm/Desktop/GDFM-CHF/aux_functions/dcc_estimation2.R")
source("/Users/ctruciosm/Desktop/GDFM-CHF/aux_functions/loglik_dcc2.R")
source("/Users/ctruciosm/Desktop/GDFM-CHF/aux_functions/R_dcc.R")
source("/Users/ctruciosm/Desktop/GDFM-CHF/aux_functions/grad_dcc2.R")

# Reading data
Y = fread("/Users/ctruciosm/Desktop/GDFM-CHF/Data/retornos_APP3.txt")[,-1]
n=dim(Y)[1]; p=dim(Y)[2]
#Setting parameters for the RW scheme.
W = 750; L = n;  WR = L-750
#Consider the two groups identified using the entire data
# 591-635  and 640-653
index = 1:p
index2 = c(c(591,635,640,653),index[-c(591,635,640,653)])
Hvector = matrix(NA,nrow=p*p,ncol=WR)
for (l in 1:50){
  print(l)
  H = matrix(0,p,p)
  h = c()
  rwY = scale(Y[l:(l+750-1),],scale = FALSE)
  Trans=segmentVOL(rwY,5)
  B = Trans$B
  Bperm = B[index2,]
  Xord = as.matrix(rwY)%*%t(Bperm)
  ##### FROM HERE #####
  # First groups: 2 elements
  GARCHaux1 = garch(Xord[,1], order = c(1,1),trace=FALSE)
  if(sum(GARCHaux1$coef[2:3])<1){
      H[1,1] = sqrt(GARCHaux1$coef[1] + GARCHaux1$coef[2]*Xord[750,1]^2+GARCHaux1$coef[3]*GARCHaux1$fitted.values[750,1]^2)
      res1 = c(as.vector(GARCHaux1$residuals[-1]),0)
  }else{
      param = GARCHestim(Xord[,1])
      h[1] = var(Xord[,1])
      for (ii in 2:751){
          h[ii] = param[1]+param[2]*Xord[ii-1,1]^2+param[3]*h[ii-1]
      }
      H[1,1] = sqrt(h[751])
      res1 = c(Xord[2:750,1]/sqrt(h[2:750]),0)
  }
  
  GARCHaux2 = garch(Xord[,2], order = c(1,1),trace=FALSE)
  if(sum(GARCHaux2$coef[2:3])<1){
      H[2,2] = sqrt(GARCHaux2$coef[1] + GARCHaux2$coef[2]*Xord[750,2]^2+GARCHaux2$coef[3]*GARCHaux2$fitted.values[750,1]^2)
      res2 = c(as.vector(GARCHaux2$residuals[-1]),0)
  }else{
      param = GARCHestim(Xord[,2])
      h[1] = var(Xord[,2])
      for (ii in 2:751){
          h[ii] = param[1]+param[2]*Xord[ii-1,2]^2+param[3]*h[ii-1]
      }
      H[2,2] = sqrt(h[751])
      res2 = c(Xord[2:750,2]/sqrt(h[2:750]),0)
  }
  
  stdres1 = cbind(res1,res2)
  parDCC1 = dcc.estimation2(stdres1[1:749,], para = c(0.05,0.9), gradient = 1)$par
  if(sum(parDCC1)>1){
      parDCC1 = dcc.estimation2(stdres1[1:749,], para = c(0.1,0.1), gradient = 1)$par
  }
  uncR <- cov(stdres1[1:749,])
  out <- .Call("dcc_est", stdres1, uncR, parDCC1[1], parDCC1[2])
  R1 = out[[1]]
  H[1:2,1:2] = H[1:2,1:2]%*%matrix(R1[750,],2,2)%*%H[1:2,1:2]
  rm(R1,out,uncR,res1,res2,stdres1,parDCC1,GARCHaux1,GARCHaux2)
  
  
  # Second groups: 2 elements
  GARCHaux3 = garch(Xord[,3], order = c(1,1),trace=FALSE)
  if(sum(GARCHaux3$coef[2:3])<1){
      H[3,3] = sqrt(GARCHaux3$coef[1] + GARCHaux3$coef[2]*Xord[750,3]^2+GARCHaux3$coef[3]*GARCHaux3$fitted.values[750,1]^2)
      res3 = c(as.vector(GARCHaux3$residuals[-1]),0)
  }else{
      param = GARCHestim(Xord[,3])
      h[1] = var(Xord[,3])
      for (ii in 2:751){
          h[ii] = param[1]+param[2]*Xord[ii-1,3]^2+param[3]*h[ii-1]
      }
      H[3,3] = sqrt(h[751])
      res3 = c(Xord[2:750,3]/sqrt(h[2:750]),0)
  }
  
  GARCHaux4 = garch(Xord[,4], order = c(1,1),trace=FALSE)
  if(sum(GARCHaux4$coef[2:3])<1){
      H[4,4] = sqrt(GARCHaux4$coef[1] + GARCHaux4$coef[2]*Xord[750,4]^2+GARCHaux4$coef[3]*GARCHaux4$fitted.values[750,1]^2)
      res4 = c(as.vector(GARCHaux4$residuals[-1]),0)
  }else{
      param = GARCHestim(Xord[,4])
      h[1] = var(Xord[,4])
      for (ii in 2:751){
          h[ii] = param[1]+param[2]*Xord[ii-1,4]^2+param[3]*h[ii-1]
      }
      H[4,4] = sqrt(h[751])
      res4 = c(Xord[2:750,4]/sqrt(h[2:750]),0)
  }
  
  stdres2 = cbind(res3,res4)
  parDCC2 = dcc.estimation2(stdres2[1:749,], para = c(0.05,0.9), gradient = 1)$par
  if(sum(parDCC2)>1){
      parDCC2 = dcc.estimation2(stdres2[1:749,], para = c(0.1,0.1), gradient = 1)$par
  }
  uncR <- cov(stdres2[1:749,])
  out <- .Call("dcc_est", stdres2, uncR, parDCC2[1], parDCC2[2])
  R2 = out[[1]]
  H[3:4,3:4] = H[3:4,3:4]%*%matrix(R2[750,],2,2)%*%H[3:4,3:4]
  rm(R2,out,uncR,res3,res4,stdres2,parDCC2,GARCHaux3,GARCHaux4)
  
  
  # Componentes GARCH
  for (i in 5:p){
      GARCHaux = garch(Xord[,i], order = c(1,1),trace=FALSE)
      if(sum(GARCHaux$coef[2:3])<1){
          H[i,i] = GARCHaux$coef[1] + GARCHaux$coef[2]*Xord[750,i]^2+GARCHaux$coef[3]*GARCHaux$fitted.values[750,1]^2
      }else{
          param = GARCHestim(Xord[,i])
          h[1] = var(Xord[,i])
          for (ii in 2:751){
              h[ii] = param[1]+param[2]*Xord[ii-1,i]^2+param[3]*h[ii-1]
          }
          H[i,i] = h[751]
      }
  }
  ##### TILL HERE #####
  # Cov(Xt/t-1) = H
  Bperminv = solve(Bperm)
  Hvector[,l] = as.vector(Bperminv%*%H%*%t(Bperminv))
}
Sigma = Hvector[,1:50]
write.table(data.frame(Sigma),"H_PCA4TS_P1_TESTE.txt",sep=",",dec=".",col.names=FALSE,row.names = FALSE)

