# Vector norm
vnorm = function(v,p=2) { 
  if ( p =="I") {
    return(max(abs(v)))
  }
  else {
    return(sum(abs(v)^p)^(1/p))
  }
}


# Solves Ax = b iteratively using the Jacobi Method
# m is the maximum number of iterations: default m = 25
# p is the p value of the matrix norm
# tol is the stopping tolerance (using the relative residual norm on the backward error)
jacobi = function(A,b,m=25,x = rep(0,n),p=2,tol=0.5*10^(-6),history=FALSE) {
  n = length(b)
  if (history) {
    hist = matrix(NA,nrow=length(b),ncol=(m+1))
    hist[,1] = x
  }
  d = diag(A)
  R = A
  R[cbind((1:n),(1:n))] = 0  # allows for r to be sparse  
  steps=0
  for (j in 1:m) {
    x = (b - R %*% x)/d 
    steps = steps+1
    if (history) {hist[,(j+1)] = as.matrix(x)}
    if (vnorm(b-A%*%x,p) <= vnorm(b,p)*tol) break 
  }
  if (history) return(list(x=x,iterations=steps,history = hist[,1:(steps+1)]))
  else return(list(x=x,steps=steps))
}