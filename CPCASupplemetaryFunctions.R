source("CPCAFunctions.R")

ProjectionPolyhedral <- function(x, A, b, initialPoint){
  # Computes the  
  
  fr <- function(z){
    sum((z - x)^2)
  }
  
  grr <- function(z){
    2*(z-x)
  }
  
  return(constrOptim(theta = initialPoint, f = fr, grad = grr, ui=A, ci=b)$par)
}

ProjectionConvex <- function(x, A, b, phi, refPoint){
  
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

CPCAFull <- function(M, X, A, b){
  
  N <- nrow(X)
  d <- ncol(X)
  
  ## Calculate the mean of data
  xMean <- colMeans(X)
  
  ## Calculate the convex principal directions 
  CPDirections <- CPCA(npcs = M, X = t(X), x0 = xMean, A = A, b = b, h = 1e-10)
  
  ## Calculate the convex principal component scores and finite-representation of each data
  CPCScores <- matrix(0, nrow = N, ncol=M)
  finiteReps <- matrix(0, nrow=N, ncol=d)
  for(i in 1:N){
    resultProjection <- ProjectionConvex(x = X[i,], A = A, b = b, phi = as.matrix(CPDirections), refPoint = xMean)
    CPCScores[i, ] <- resultProjection$coefs
    finiteReps[i, ] <- resultProjection$projs
  }
  
  ## Calculate the variance for each convex principal direction
  varianceVec <- rep(0, M)
  for(j in 1:M){
    varianceVec[j] <- apply(CPCScores, MARGIN = 2, var)
  }
  
  resultCPCA <- list()
  resultCPCA$mean <- xMean
  resultCPCA$directions <- CPDirections
  resultCPCA$scores <-CPCScores
  resultCPCA$variance <- varianceVec
  resultCPCA$finiteReps <- finiteReps
  return(resultCPCA)
}



ConvexCVE <- function(maxComp, X, A, b){
  
  N <- nrow(X)
  d <- ncol(X)
  
  # Calculate the total variance
  meanX <- colMeans(X)
  TV <- mean(apply(X - meanX, MARGIN = 1, function(x){sum(x^2)}))
  
  # Calculate the cumulative variance explained
  CVE <- rep(0, maxComp)
  for(M in 1:maxComp){
    resultCPCA <- CPCAFull(M, X, A, b)
    finiteReps <- resultCPCA$finiteReps
    variationByM <- mean(apply(finiteReps - meanX, MARGIN = 1, function(x){sum(x^2)}))
    CVE[M] <- variationByM / TV
  }
  
  return(CVE)
}






