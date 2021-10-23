uncons.optim.mvp = function(H){
  ones = rep(1,dim(H)[1])
  return(as.numeric(solve(H)%*%ones)/(t(ones)%*%solve(H)%*%ones))
}


func.mvp = function(omega,H){
  #omega = par# c(par,1-sum(par))
 # if (min(omega)>=0 && max(omega)<= 1 && sum(omega)==1){
    Func = t(omega)%*%H%*%omega
  #} else{
  #  K = dim(H)[1]
  #  Func = t(rep(1/K,K))%*%H%*%rep(1/K,K)
  #}
  return(Func)
}
equal <- function(omega, H) {
  sum(omega)
}

optim.mvp = function(H){
  K = dim(H)[1]
  pesosini = rep(1/K,K)
  paraux2 = solnp(pars = pesosini, fun = func.mvp, eqfun=equal,eqB=1,LB=rep(0,K), UB=rep(1,K), H = H)$pars
  #paraux2 = solnp(pars = pesosini, fun = func.mvp, eqfun=equal,eqB=1,LB=rep(-Inf,K), UB=rep(Inf,K), H = H)$pars
  return (paraux2) 
}




