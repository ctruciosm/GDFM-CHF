NumComp = function(evalues){
  N = length(evalues)
  R = matrix(NA,ncol=1,nrow=2)
  if(is.null(evalues)){
    stop("Vector is null")
  }
  if(is.matrix(evalues)){
    stop("Input is a matrix, vector is required")
  }
  R[1,1] = sum(evalues>mean(evalues))
  R[2,1] = which(evalues[1:(N-1)]/evalues[2:N]==max(evalues[1:(N-1)]/evalues[2:N]))
  return(R)
}