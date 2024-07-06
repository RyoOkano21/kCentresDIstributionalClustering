DyadicPartition <- function(a,b,n){
  # Generates the nth Dyadic partition of the interval [a,b]
  
  D <- a+ ((b-a) * ((1/2)^n)) * seq(0,2^n,1)
  
  return(D)
}


ConstraintMatrixGPCA <- function(n){
  # Generates the constraint matrix for the (discrete) GPCA problem
  # defined on the nth dyadic partition of the interval [a,b]
  
  A <- matrix(0, nrow=2^n+1,ncol=2^n)
  for(i in seq(1,2^n+1,1)){
    if(i > 1){
      if(i<2^n+1){
        A[i,i] <- 1
        A[i,i-1]<--1
      }else{
        A[i,i-1] <- -1
      }
    }else{
      A[i,i] <- 1
    }
  }
  return(A)
}


ConstraintVectorGPCA <-function(a,b,n){
  # Generates the constraint vector for the (discrete) GPCA problem
  # defined on the nth dyadic partition of the interval [a,b]
  c <- DyadicPartition(a, b, n)[-1]
  v <- rep(0, 2^n+1)
  v[1] <- a-c[1]
  for(i in 2:(2^n)){
    v[i] <- c[i-1] - c[i]
  }
  v[2^n+1] <- c[2^n]-b
  return(v)
}