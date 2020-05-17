Cond = function(A,p=2) {
  if (p == 2) {  # by default use the 
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

