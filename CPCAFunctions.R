library(nloptr)
library(Matrix)
library(profvis)
library(Rcpp)
library(RcppParallel)
library(foreach)
library(doParallel)
library(pracma)
library(parallel)
library(readr)
library(dplyr)
library(rbenchmark)
library(stats)
library(MASS)
library(fossil)

sourceCpp("HelperFunctions.cpp")


eval_CPCA_obj_closure <- function(theta,B,x0,A,b,X) {
  # Creates a function of theta for optimization
  function(theta) eval_CPCA_obj(theta,B,x0,A,b,X)
}

eval_grad_CPCA_obj_closure <- function(theta,h,B,x0,A,b,X) {
  # Creates a function of theta for optimization
  function(theta) eval_grad_CPCA_obj(theta,h,B,x0,A,b,X)
}

CPdirection <- function(npcs,X,x0,A,b,h){
  #Computes the convex principal directions
  
  # npcs: number of principal components requested
  # X: data matrix (each row corresponds to a variable, and column corresponds to a case)
  # x0: reference element
  # A: constraint matrix
  # b: constraint vector
  # h: central difference step size
  
  # Extract Dimension
  d<-nrow(X)
  
  # Allocate Matrix
  PCs<-matrix(0,ncol=npcs,nrow=d)
  
  # Optimization specifications
  opts <- list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8)
  
  for(i in seq(1,npcs,1)){
    # Get principle component estimate, projected data, and Basis Matrix
    if(i>1){
      X1<-project_orthogonal_complement(X-x0,as.matrix(PCs[,1:(i-1)]))
      B<-pracma::nullspace(t(PCs[,1:(i-1)]))
      pca <- prcomp(t(X1),center=FALSE,scale=FALSE)
      p<-pca$rotation[,1]
      p0<-ls_fit(B,p)
    }else{
      B<-eye(d)
      pca <- prcomp(t(X-x0),center=FALSE,scale=FALSE)
      p<-pca$rotation[,1]
      p0<-p
    }
    
    # Get implied theta
    theta0 <- as.vector(spherical_coords(p0))
    
    # Get functions for optimization routine
    eval_f<-eval_CPCA_obj_closure(theta,B,x0,A,b,X)
    eval_grad_f <- eval_grad_CPCA_obj_closure(theta,h,B,x0,A,b,X)
    
    status <- paste0("Optimization Running for PC: ", i, " ...")
    
    # Run optimization
    res <- nloptr(
      x0          = theta0,
      eval_f      = eval_f,
      eval_grad_f = eval_grad_f,
      opts        = opts
    )
    theta_sol<-res$solution
    status <- paste0("Optimization for PC: ", i, " Complete.")
    
    # Store solution
    PCs[,i]<-B%*%undo_spherical_coords(theta_sol)
  }
  
  return(PCs)
}


ProjectionPolyhedral <- function(x, A, b, initialPoint){
  # Computes the projection of x onto the polyhedral domain defined by A and b
  
  fr <- function(z){
    sum((z - x)^2)
  }
  
  grr <- function(z){
    2*(z-x)
  }
  
  return(constrOptim(theta = initialPoint, f = fr, grad = grr, ui=A, ci=b)$par)
}

ProjectionConvex <- function(x, A, b, phi, refPoint){
  # Computes the projection of x onto the intersection of the polyhedral domain 
  # and space spaned by the vectors in phi
  
  M <- ncol(phi)
  
  # Calculate the projection of x onto the space spaned by the vectors in phi
  coefs <- as.vector(t(x-refPoint) %*% phi)
  xProj <- refPoint + phi %*% coefs
  
  # If the projection is not contained in the convex set, do further projection
  if(!(all(A %*% xProj >= b))){
    constraintMat <- A %*% phi
    constraintVec <- b - A %*% refPoint
    coefs <- ProjectionPolyhedral(x = coefs, A = constraintMat, b = constraintVec-10^(-10), initialPoint = rep(0, M))
    xProj <- refPoint + phi %*% coefs
  }
  
  resultProjection <- list()
  resultProjection$coefs <- coefs
  resultProjection$projs <- xProj
  return(resultProjection)
}

ConvexPCA <- function(M, X, A, b){
  # Performs convex PCA 
  
  # Input 
  # M: number of principal components requested
  # X: data matrix (each row corresponds to a case, and column corresponds to a row)
  # A: constraint matrix
  # b: constraint vector
  
  # Output
  # mean: the mean vector 
  # directions: a matrix containing the convex principal directions 
  # scores: a matrix containing the convex principal component scores 
  # variances: the variance for each 
  # finiteReps: the finite-dimensional representations 
  
  N <- nrow(X)
  d <- ncol(X)
  
  # Calculate the mean of data
  xMean <- colMeans(X)
  
  # Calculate the convex principal directions 
  CPDirections <- CPdirection(npcs = M, X = t(X), x0 = xMean, A = A, b = b, h = 1e-10)
  
  # Calculate the convex principal component scores and finite-representation of each data
  CPCScores <- matrix(0, nrow = N, ncol=M)
  finiteReps <- matrix(0, nrow=N, ncol=d)
  for(i in 1:N){
    resultProjection <- ProjectionConvex(x = X[i,], A = A, b = b, phi = as.matrix(CPDirections), refPoint = xMean)
    CPCScores[i, ] <- resultProjection$coefs
    finiteReps[i, ] <- resultProjection$projs
  }
  
  ## Calculate the variance for each convex principal direction
  varianceVec <- rep(0, M)
    varianceVec <- apply(CPCScores, MARGIN = 2, var)
  
  resultCPCA <- list()
  resultCPCA$mean <- xMean
  resultCPCA$directions <- CPDirections
  resultCPCA$scores <-CPCScores
  resultCPCA$variances <- varianceVec
  resultCPCA$finiteReps <- finiteReps
  return(resultCPCA)
}



ConvexCVE <- function(maxComp, X, A, b){
  # computes the cumulative variance explained
  
  N <- nrow(X)
  d <- ncol(X)
  
  # Calculate the total variance
  meanX <- colMeans(X)
  varVec <- rep(0, N)
  for(i in 1:N){
    varVec[i] <- sum((X[i, ] - meanX)^2)
  }
  TV <- mean(varVec)
  
  # Calculate the cumulative variance explained
  CVE <- rep(0, maxComp)
  for(M in 1:maxComp){
    resultCPCA <- ConvexPCA(M, X, A, b)
    finiteReps <- resultCPCA$finiteReps
    varVec <- rep(0, N)
    for(i in 1:N){
      varVec[i] <- sum((finiteReps[i, ] - meanX)^2)
    }
    variationByM <- mean(varVec)
    CVE[M] <- variationByM / TV
  }
  
  return(CVE)
}





