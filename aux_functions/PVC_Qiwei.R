comVolQ = function (rtn, m = 10) {
  x1 = scale(rtn, center = T, scale = F)
  nT = dim(x1)[1]
  k = dim(x1)[2]
  CovData = t(x1)%*%x1/nT
  A = diag(x1%*%t(x1))
  Sum1 = matrix(0, k, k)
  for (kk in 1:m) {
    Sum2 = matrix(0, k, k)
    for (tt in 1:nT){
      sign = (A<= A[tt]) 
      signnew = c(rep(0,kk),sign[1:(nT-kk)])
      temp1 = x1*matrix(signnew,ncol=k,nrow=nT)
      temp2 = (t(temp1)%*%temp1-sum(signnew)*CovData)/(nT-kk)
      Sum2 = Sum2 + temp2%*%temp2/nT
    }
  Sum1 = Sum1 + Sum2
  }
  m2 = eigen(Sum1)   
  Valu = m2$values
  Prop = Valu/sum(Valu)
  Vec = m2$vectors
  comVolQ <- list(residuals = x1, values = m2$values, vectors = m2$vectors,COV = CovData)
}