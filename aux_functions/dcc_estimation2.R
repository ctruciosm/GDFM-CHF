#*****************************************************************************************************************
# The 2nd step DCC estimation.
dcc.estimation2 <- function(dvar, para, gradient=0){ # dvar must be standardised residuals
   resta <- rbind(c(-1, -1), diag(2))
   restb <- c(-0.9999, 0.0001, 0.0001)
   if(gradient!=0){
      step2 <- constrOptim(theta=para, f=loglik.dcc2, grad=grad.dcc2, ui=resta, ci=restb, mu=1e-5, dvar=dvar)
   } else {
      step2 <- constrOptim(theta=para, f=loglik.dcc2, grad=NULL, ui=resta, ci=restb, mu=1e-5, dvar=dvar)
   }
   step2
}



loglik.garch <- function(param, dvar){
  nobs <- length(dvar)
  h = c()
  h[1] = var(dvar)
  for(i in 2:nobs){
    h[i] = param[1] + param[2]*dvar[i-1]^2+param[3]*h[i-1]
  }
  z <- dvar/sqrt(h)
  lf <- -0.5*sum(log(h))-0.5*sum(z^2)
    return(-lf)
}



 GARCHestim = function(datos){
   resta <- matrix(c(1,0,0,0,1,0,0,0,1,0,-1,-1),ncol=3,byrow=TRUE)
   restb <- c(0.0001, 0.0001, 0.0001,-0.9999)
   step1 <- constrOptim(theta=c(0.001,0.05,0.8), f=loglik.garch, grad=NULL,ui=resta, ci=restb, mu=1e-5, dvar=datos)$par
   return(step1)
 }

### Rascunhos
 
 # gridgarch = function(dvar){
 #   param = c(0.001,0.05,0.8)
 #   Fun = loglik.garch(param,dvar)
 #   omega = seq(0.001,0.2,length.out=5)
 #   alpha = seq(0.001,0.5,length.out=7)
 #   beta = seq(0.001,0.99,length.out=10)
 # for (i in omega){
 #   for (j in alpha){
 #     for (k in beta){
 #       if(j+k<1){
 #         Fun2 = loglik.garch(c(i,j,k),dvar)
 #         if(Fun2<Fun){
 #           param= c(i,j,k)
 #           Fun=Fun2
 #         }
 #       }
 #     }
 #   }
 # }
 #   return(param)
 # }
 
 # Alternative optimization based on ccgarch2 package
 
 # loglik2.dcc <- function(param, data){        # data is the standardised residuals
 #   if(is.zoo(data)) data <- as.matrix(data)
 #   nobs <- dim(data)[1]
 #   ndim <- dim(data)[2]
 #   DCC <- dcc.est(data, param)$DCC
 #   #lf1 <- dcc.ll2(DCC, data)    
 #      lf1 <- numeric(nobs)
 #      for( i in 1:nobs){
 #        R1 <- matrix(DCC[i,], ndim, ndim)
 #        invR1 <- solve(R1)
 #        lf1[i] <- 0.5*(log(det(R1)) + sum(data[i,]*crossprod(invR1, data[i,])))  # negative of the log-likelihood function
 #      }
 #   sum(lf1)
 # }
 # 
 # 
 # inEQdcc2 <- function(param, data){
 #   param[1] + param[2]
 # }
 # 
 # dcc.2stg <- function(data, para){ # data must be standardised residuals
 #   LB <- rep(0.0001, 2)     # lower bound of the parameter
 #   UB <- rep(0.9999, 2)  
 #   ineqLB <- 0.0001         # lower bound of the stationarity
 #   ineqUB <- 0.9999         # upper bound of the stationarity
 #   
 #   suppressWarnings(
 #     step2 <- solnp(pars=para, fun = loglik2.dcc,
 #                    ineqfun = inEQdcc2, ineqLB = ineqLB, ineqUB = ineqUB,
 #                    LB = LB, UB = UB,
 #                    control = list(trace=0, tol = 1e-6, delta = 1e-10),
 #                    data = data))
 #   step2
 # }
 
 
 
 # First step
