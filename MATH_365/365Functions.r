#########################################################################
# Norms, condition numbers, and special matrices

# Vector norm
vnorm = function(v,p=2) { 
  if ( p =="I") {
    return(max(abs(v)))
  }
  else {
    return(sum(abs(v)^p)^(1/p))
  }
}

# Condition number of a matrix
Cond = function(A,p=2) {
  if (p == 2) {  # by default use the 2-norm
    s = svd(A)$d
    s = s[s>0]
    return(max(s)/min(s))
  }
  if (p == 1) {  # use the 1 norm
    Ainv = solve(A)
    return(max(colSums(abs(A)))*max(colSums(abs(Ainv))))
  }
  if (p == 'I') {   # use the infinity norm
    Ainv = solve(A)
    return(max(rowSums(abs(A)))*max(rowSums(abs(Ainv))))
  }
}

# Create a tridiagonal matrix. The diagonal is d, the first superdiagonal is u and the first subdiagonal is l
TriDiag = function(d,l,u, sparse=TRUE) {
  n = length(d)
  if ( length(l) != (n-1) || length(u) != (n-1) ) 
    stop('vectors not of correct tridiagonal length')
  if (sparse) {
    D = Matrix(0,nrow=n,ncol=n,sparse=TRUE)
    diag(D) = d
  }
  else  D = diag(d)
  D[cbind((1:(n-1)),(2:n))] = u
  D[cbind((2:n),(1:(n-1)))] = l
  return(D)
}


#########################################################################
# Chapter 0 and 1: Function evaluation and root-finding methods

# Horner's method performs polynomial evaluation by nesting
# It takes as input a value of x and the coefficients of a polynomial 
# ordered from the constant term up to the nth term
# It returns p(x)

Horner = function(coeffs,x,b=rep(0,length(coeffs))) {
  y = coeffs[length(coeffs)]
  for (i in (length(coeffs)-1):1) {
    y = y*(x-b[i])+coeffs[i] 
  }
  return(y)
}

# Wilkinson's polynomial
Wilkinson = function(x) {
  (-20 + x)*(-19 + x)*(-18 + x)*(-17 + x)*(-16 + x)*(-15 + x)*(-14 + x)*(-13 + x)*
    (-12 + x)*(-11 + x)*(-10 + x)*(-9 + x)*(-8 + x)*(-7 + x)*(-6 + x)*(-5 + x)*(-4 + x)*
    (-3 + x)*(-2 + x)*(-1 + x)
}

W = function(x){
  2432902008176640000 - 8752948036761600000*x + 13803759753640704000*x^2 - 
    12870931245150988800*x^3 + 8037811822645051776*x^4 - 3599979517947607200*x^5  + 
    1206647803780373360*x^6 - 311333643161390640*x^7 + 63030812099294896*x^8 - 
    10142299865511450*x^9 + 1307535010540395*x^10 - 135585182899530*x^11 + 
    11310276995381*x^12 - 756111184500*x^13 + 40171771630*x^14 -1672280820*x^15 + 
    53327946*x^16 - 1256850*x^17 + 20615*x^18 - 210*x^19 + x^20
}


# Finite difference derivative
FDDeriv = function(f,delta=.000001){ 
  function(x){(f(x+delta) - f(x-delta))/(2*delta)}
}

# Bisection method

bisect = function(f,interval,tol=0.5*10^-10,maxiters=100){
  # pre-allocate memory
  left.brackets = rep(NA, maxiters)
  right.brackets = rep(NA, maxiters)
  approximations = rep(NA, maxiters)
  fleft = rep(NA, maxiters)
  fright = rep(NA, maxiters)
  
  left.brackets[1]=interval[1]
  right.brackets[1]=interval[2]
  fleft[1]=f(left.brackets[1])
  fright[1]=f(right.brackets[1])
  j = 1                                
  if (sign(fleft[1])*sign(fright[1]) >= 0 ) 
    stop("f(a)f(b)<0 not satisfied")
  while ((right.brackets[j]-left.brackets[j])/2 > tol & j <= maxiters ) {
    approximations[j] = (left.brackets[j] + right.brackets[j])/2
    fc=f(approximations[j])
    if (fc == 0)  break
    if (sign(fc)*sign(fleft[j]) < 0) {
      right.brackets[j+1] = approximations[j]
      fright[j+1] = fc
      left.brackets[j+1]=left.brackets[j]
      fleft[j+1] = fleft[j]
    }
    else {
      right.brackets[j+1] = right.brackets[j] 
      fright[j+1] = fright[j]
      left.brackets[j+1]=approximations[j]
      fleft[j+1] = fc
    }
    j = j+1
  }
  return(list(root = (left.brackets[j]+right.brackets[j])/2, approximations = approximations[!is.na(approximations)],left.brackets = left.brackets[!is.na(left.brackets)],right.brackets = right.brackets[!is.na(right.brackets)],fleft = fleft[!is.na(fleft)],fright = fright[!is.na(fright)])) 
}

# a shorter version of bisect that just returns the root
bisect.short = function(f,interval,tol=0.5*10^-10,verbose=FALSE){
  a = interval[1]
  b = interval[2]
  fa = f(a)
  fb = f(b)
  j = 0                                # counter for verbose printing
  if (sign(fa)*sign(fb) >= 0 ) 
    stop("f(a)f(b)<0 not satisfied")
  while ((b-a)/2 > tol) {
    c = (a + b)/2
    fc = f(c)
    j = j+1
    if (verbose==TRUE) {               # show values at each step
      print(j) 
      print(c(a,c,b,b-a))
    }
    if (fc == 0)  break
    if (sign(fc)*sign(fa) < 0) {
      b = c
      fb = fc
    }
    else {
      a = c
      fa = fc
    }
  }
  return((a+b)/2) 
}

# Secant method
secant = function(f, a, b, tol = 1e-05, maxiters = 12) {
  history = rep(NA, maxiters)
  history[1] = a
  history[2] = b
  for (k in 3:maxiters) {
    x1 = history[k-2]
    x2 = history[k-1]
    newx = x2 - (f(x2)*(x2-x1))/(f(x2)-f(x1))
    history[k] = newx
    change = newx - history[k-1]
    if (abs(change) < tol * max(abs(history[k - 1]), tol)) break
  }
  return(list(root = newx, history = history[!is.na(history)]))
}

# Inverse quadratic interpolation method
iqi = function(f, a, b,c, tol = 1e-05, maxiters = 15) {
  history = rep(NA, maxiters)
  history[1] = a
  history[2] = b
  history[3] = c
  for (k in 4:maxiters) {
    x0 = history[k-3]
    x1 = history[k-2]
    x2 = history[k-1]
    q1=f(x0)/f(x1)
    q2=f(x2)/f(x1)
    q3=f(x2)/f(x0)
    newx = x2-(q2*(q2-q1)*(x2-x1)+(1-q2)*q3*(x2-x0))/((q1-1)*(q2-1)*(q3-1))
    history[k] = newx
    change = newx - history[k-1]
    if (abs(change) < tol * max(abs(history[k - 1]), tol)) break
  }
  return(list(root = newx, history = history[!is.na(history)]))
}

#########################################################################
# Chapter 2: Methods to solve Ax=b

# Part I: Direct Methods
########################

# Naive Gaussian elimination without pivoting
eliminate = function(A,tol=10^-8) {
  n = nrow(A)
  for ( j in 1:(n-1) ) {
    pivot = A[j,j]
    if (abs(pivot) < tol) stop('zero pivot encountered')
    for ( i in (j+1):n ) {
      A[i,] = A[i,] - A[i,j]/pivot * A[j,]
    }
  }
  return(A)
}

MySolve = function(A,b,tol=10^-8) {
  n = nrow(A)
  A = cbind(A,b)     # append b to the right of A
  A = eliminate(A)  # Gaussian elimination below the main diagonal
  b = A[,n+1]        # Get the row reduced version of b back
  x = rep(0,n)       # pre-allocate a vector x with 0s in it.
  x[n] = b[n]/A[n,n] # Fill in the nth value of A
  for (j in (n-1):1 ) {
    back = 0                  #  compute the back substitution part with a for loop
    for (k in (j+1):n ) {
      back = back + A[j,k]*x[k]
    }
    x[j] = (b[j] - back)/A[j,j]
  }
  return(x)
}

# The following version is vectorized and does not
# require an inner for loop to compute the back substitution

MySolveV = function(A,b,tol=10^-8) {
  n = nrow(A)
  A = cbind(A,b)     # append b to the right of A
  A = eliminate(A)  # Gaussian elimination below the main diagonal
  b = A[,n+1]        # Get the row reduced version of b back
  x = rep(0,n)       # pre-allocate a vector x with 0s in it.
  x[n] = b[n]/A[n,n] # Fill in the nth value of A
  for (j in (n-1):1 ) {
    back = A[j,((j+1):n)] %*% x[(j+1):n]   # compute the back sub part with a dot product
    x[j] = (b[j] - back)/A[j,j]
  }
  return(x)
}

# LU without pivoting
myLU = function(A,tol=10^-8) {
  n = nrow(A)
  L = diag(x=1,nrow=n)  # start with L = identity matrix
  U=A
  for ( k in 1:(n-1) ) {
    pivot = U[k,k]
    if (abs(pivot) < tol) stop('zero pivot encountered')
    for ( j in (k+1):n ) {
      mult =  U[j,k]/pivot
      U[j,] = U[j,] - mult * U[k,]
      L[j,k] = mult   #  store the multiplier in LA
    }
  }
  return(list(L=L,U=U))
}

myLUFast = function(A,tol=10^-8) {
  n = nrow(A)
  L = diag(x=1,nrow=n)  # start with L = identity matrix
  U=A
  for ( k in 1:(n-1) ) {
    pivot = U[k,k]
    if (abs(pivot) < tol) stop('zero pivot encountered')
    mults=U[(k+1):n,k]/pivot
    U[(k+1):n,k]=0
    U[(k+1):n,(k+1):n]=U[(k+1):n,(k+1):n]-mults%o%U[k,(k+1):n]
    L[(k+1):n,k]=mults
  }
  return(list(L=L,U=U))
}

myLUSolve=function(L,U,b,tol=1e-10){
  
  n=nrow(L)
  
  # First solve Ly=b 
  y = rep(0,n)        # pre-allocate a vector y with 0s in it.
  if(abs(L[1,1])<tol) stop('There is a zero on the diagonal of L')
  y[1] = b[1]/L[1,1]  # Fill in the 1st value of y
  for (j in 2:n ) {
    if(abs(L[j,j])<tol) stop('There is a zero on the diagonal of L')
    y[j] = (b[j] - L[j,1:(j-1)]%*%y[1:(j-1)])/L[j,j]
  }
  
  # Then solve Ux=y 
  x = rep(0,n)        # pre-allocate a vector x with 0s in it.
  if(abs(U[n,n])<tol) stop('There is a zero on the diagonal of U')
  x[n] = y[n]/U[n,n]  # Fill in the nth value of x
  for (j in (n-1):1 ) {
    if(abs(U[j,j])<tol) stop('There is a zero on the diagonal of U')
    x[j] = (y[j] - U[j,(j+1):n]%*%x[(j+1):n])/U[j,j]
  }
  return(x)
}

# Part II: Iterative Methods
############################

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

# Gauss-Seidel Method
GaussSeidel = function(A,b,k=10) {
  n = nrow(A)
  x = rep(0,n)
  d = diag(A)
  L = A
  U = A
  U[lower.tri(A,diag=TRUE)] = 0  # U is the upper triangular part of A
  L[upper.tri(A,diag=TRUE)] = 0  # L is the lower triangular part of A
  for (j in 1:k) {
    bu = b - U %*% x      
    for (i in 1:n)  
      x[i] = (bu[i] - L[i,] %*% x)/d[i]    # successivly update x as we use it.
  }
  return(x)
}

# Conjugate Gradient Method
# Inputs: symm pos def matrix A, rhs b, number of steps n, 
# x is the initial guess, 
# tol is a stopping condition on the size of the residual
# history gives a full history
# Output: solution x to Ax=b 
ConjGrad = function(A,b,x = rep(0,length(b)),tol=1e-18,m = length(b),history=FALSE) {
  n = length(b)
  r = b - A %*% x
  d = r
  if (history) {
    hist = matrix(NA,nrow=n,ncol=m+1)
    hist[,1] = x
  }
  for (i in 1:m ) {
    if (max(abs(r)) < tol) break
    steps = i
    alpha = (t(r) %*% r)/(t(d) %*% A %*% d)  # step length
    x = x + alpha[1,1] * d                  # take step 
    if (history) {hist[,i+1] = x[,1]}
    rold = r
    r = rold - alpha[1,1] * A %*% d         # new residual
    beta = (t(r) %*% r)/(t(rold) %*% rold)   # improvement this step
    d = r + beta[1,1]*d
  }
  if (history) {return(list(x=x,history=hist))}
  else {return(list(x=x,steps=steps))}
}


# Unit Circle Mapping Under a Matrix A
UnitCircleMap = function(A, p = 2) {
  
  if (p == 2) {
    t = seq(0, 2 * pi, len = 1000)
    x = cos(t)
    y = sin(t)
    nn = base::norm(A, type = "2")
    ni = as.character(p)
  } else if (p == 1) {
    # Insert some code here
    nn = base:norm(A, type = "1")
    ni = as.character(p)
  } else if (p == "I") {
    # Insert some code here
    nn = base::norm(A, type = "I")
    ni = "Inf"
  }
  
  pts = A %*% t(cbind(x, y))
  
  newx = pts[1, ]
  newy = pts[2, ]
  
  M = max(c(newx, newy, 1.5))
  m = min(c(newx, newy, -1.5))
  
  plot(x, y, type = "l", col = "black", xlim = c(m, M), ylim = c(m, M), xlab = "x1", 
       ylab = "x2")
  lines(newx, newy, col = "red")
  
  normSize = sprintf("||A||_%s = %1.2f", ni, nn)
  text(-0.7 * M, 0.9 * M, normSize)
}



#########################################################################
# Chapter 3: Interpolation and Function Approximation

# Produces the Vandermonde matrix
Vandermonde = function(x) {
  n = length(x)
  
  # Method 1: a for loop
  #  V = matrix(0,nrow=n, ncol=n)
  #  for (j in 1:n)
  #    V[,j] = x^(j-1)
  
  # Method 2: an outer product (a little cleaner)
  V=outer( x, 0:(n-1), `^` )
  return(V)
}

# Performs Vandermonde interpolation
VandermondeInterpolator = function(x,y) {
  n = length(x)
  V = Vandermonde(x)
  c = solve(V,y)
  p = function(z) { 
    Horner(c,z)
    #c %*% z^(0:(length(c)-1)) this alternative is just fine as well
  }  
  xx = seq(min(x)-1,max(x)+1,length=1000)
  yy=p(xx)
  plot(xx,yy,type='l',col='blue',xlab="x",ylab="y")
  points(x,y,pch=20,col='red')  
  return(list(coeffs=c,interpolant=p))
}


NewtonDD = function(x,y,ddtable=FALSE) {
  n = length(x)
  v = matrix(0,nrow=n,ncol=n)
  v[,1] = y
  for (i in 2:n) {
    for (j in 1:(n+1-i)) {
      v[j,i] = (v[j+1,i-1]-v[j,i-1])/(x[j+i-1]-x[j])
    }
  }
  if (ddtable) print(v)
  return(v[1,])
}

# This function takes in some interpolation points (the coordinates of which are xi and yi) and creates
# a polynomial p(x) passing through these points. The third input xx is a sequence of points at which
# this interpolating polynomial is evaluated. The output is the series of points y=p(xx); i.e. the value 
# of the polynomial at each of the points in x. 

Interpolator = function(xi,yi,xx,Itype = 'NewtonDD') {
  if (Itype == 'NewtonDD') {
    c = NewtonDD(xi,yi)
    return(Horner(c,xx,xi))
  }
  if (Itype == 'Vandermonde') {
    V = Vandermonde(xi)
    c = solve(V,yi)
    if (length(xx)==1)
      return (c %*% xx^(0:(length(c)-1)))
    else
      return(Horner(c,xx))
  }
}


# Bezier Curves

# This function returns points along a Bezier spline with given
# endpoints (x1,y1) and (x4,y4) and control points (x1,y2) and (x3,y3)
one.bezier = function(x1, x2, x3, x4, y1, y2, y3, y4, npts = 10) {
  bx = 3 * (x2 - x1)
  cx = 3 * (x3 - x2) - bx
  dx = x4 - x1 - bx - cx
  by = 3 * (y2 - y1)
  cy = 3 * (y3 - y2) - by
  dy = y4 - y1 - by - cy
  t = seq(0, 1, length = npts)
  x = x1 + t * (bx + t * (cx + dx * t))
  y = y1 + t * (by + t * (cy + dy * t))
  return(list(x = x, y = y))
}

# This function draws a single Bezier curve with given
# endpoints (x1,y1) and (x4,y4) and control points (x2,y2) and (x3,y3)
draw.one.bezier = function(x1, x2, x3, x4, y1, y2, y3, y4, npts = 10,new = TRUE, ...) {
  pts = one.bezier(x1, x2, x3, x4, y1, y2, y3, y4, npts = npts)
  if (new) 
    plot(pts, type = "n", 
         xlim = c(min(pts$x, x1, x2, x3, x4),max(pts$x, x1, x2, x3, x4)), 
         ylim = c(min(pts$y, y1, y2, y3, y4),max(pts$y, y1, y2, y3, y4)),...)
  lines(pts, ...)
  points(c(x1, x4), c(y1, y4), pch = 20, ...)
  points(c(x2, x3), c(y2, y3), pch = 10, ...)
  lines(c(x1, x2), c(y1, y2), lty = 3)
  lines(c(x4, x3), c(y4, y3), lty = 3)
  invisible(pts)
}

# This function draws multiple Bezier curves on the same plot given vectors of x and y points
draw.beziers = function(x1, x2, x3, x4, y1, y2, y3, y4, npts = 10,new = TRUE, ...) {
  allx = c(x1, x2, x3, x4)
  ally = c(y1, y2, y3, y4)
  xlim = c(min(allx), max(allx))
  ylim = c(min(ally), max(ally))
  if (new) 
    plot(-1:2, xlim = xlim, ylim = ylim,...,col='white')
  xpts = c()
  ypts = c()
  for (k in 1:length(x1)) {
    pts = draw.one.bezier(x1[k], x2[k], x3[k], x4[k], y1[k], 
                          y2[k], y3[k], y4[k], npts=npts,new = FALSE, xlim = xlim, ylim = ylim, 
                          ...)
    xpts = c(xpts, pts$x)
    ypts = c(ypts, pts$y)
  }
  invisible(list(x = xpts, y = ypts))
}

# This function allows the user to manipulate a Bezier curve.  You must include the manipulate package
bezier.interact = function(x1=0,y1=0,x2=10,y2=10){
  ang1 = 0;
  ang2 = 0;
  len1 = 1;
  len2 = 1;
  do.it = function(len1,len2,ang1,ang2) {
    cx1 = x1 + len1*cos(ang1)
    cy1 = y1 + len1*sin(ang1)
    cx2 = x2 + len2*cos(ang2)
    cy2 = y2 + len2*sin(ang2)
    #    foo = make.all.bezier(cbind(x1,y1,cx1,cy1,cx2,cy2,x2,y2))
    draw.one.bezier(x1,cx1,cx2,x2,y1,cy1,cy2,y2,npts=100,xlab="",ylab="")
  }
  manipulate( do.it(len1,len2,ang1*pi/180, ang2*pi/180), 
              len1=slider(0,10,step=.1,initial=1),
              ang1 = slider(-180,180,step=1,initial=-90),
              len2 = slider(0,10,step=.1, initial=1),
              ang2 = slider(-180,180,step=1,initial=90) )
}


# This function takes a matrix T of points in 
# x1, y1, x2, y2, x3, y3, x4, y4 forms and draws the Bezier curve
draw.curve = function(T, npts = 100, ...) {
  draw.beziers(T[, 1], T[, 3], T[, 5], T[, 7], T[, 2],T[, 4], T[, 6], T[, 8],npts=npts,xlab="",ylab="")
}
