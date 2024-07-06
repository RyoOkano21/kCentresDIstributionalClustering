source("CPCAFunctions.R")
source("GPCAFunctions.R")


GetDistances <- function(dataMat, listOfCPCAobjs, A, b){
  # Computes the distances between data and their projections onto the spaces spanned by
  # CPCA structures (means and convex principal directions)
  
  N <- nrow(dataMat)
  K <- length(listOfCPCAobjs)
  distanceMat <- matrix(0, nrow = N, ncol=K)
  
  for(i in 1:N){
    for(c in 1:K){
      # choose a data and CPCA structure
      x <- dataMat[i, ]
      cpcaObj <- listOfCPCAobjs[[c]]
      
      # calculate the projection of x
      refPoint <- cpcaObj$mean
      phi <- cpcaObj$direction
      xProj <- ProjectionConvex(x, A, b, phi, refPoint)$projs
      
      # calculate the distance between x and its projection
      distanceMat[i, c] <- sqrt(sum((x-xProj)^2))
    }
  }
  
  return(distanceMat)
}

DistributionCVE <- function(maxComp, ObsList, a, b, nd){
  # computes the 
  
  N <- length(ObsList)
  
  # get the Dyadic partitions of the intervals [a, b] and [0, 1]
  omegaGrid <- DyadicPartition(a, b, nd)[-1]
  pVals <- punif(q = omegaGrid, min = a, max = b)
  
  # calculate functions in the tangent space evaluated at the partition of [a, b]
  # the uniform distribution on [a, b] is used as the reference measure
  LFs <- matrix(0, nrow = N, ncol=2^(nd))
  for(i in 1:N){
    LFs[i, ] <- quantile(x = ObsList[[i]], pVals) - omegaGrid
  }
  
  # get the constraint matrix and vector
  A <- ConstraintMatrixGPCA(nd)
  b <- ConstraintVectorGPCA(a, b, nd)
  
  # get the cumulative variance explained
  resultDCVE <- ConvexCVE(maxComp, X = LFs, A, b)
  
  return(resultDCVE)
}


kCDC <- function(ObsList, a, b, nd, K, initialM, iterM, nStart, iterMax){
  # Performs the proposed k centres distributional clustering
  
  # input
  # ObsList: 
  # a: 
  # b: 
  # nd: 
  # K:
  # initialM: 
  # iterM: 
  # nStart: 
  # iterMax
  
  # output
  # cluster: 
  # cpcaList: 
  # iterToConv: 
  # prevConf: 
  # clustConf0: 
  
  N <- length(ObsList)
  
  # get the Dyadic partitions of the intervals [a, b] and [0, 1]
  omegaGrid <- DyadicPartition(a, b, nd)[-1]
  pVals <- punif(q = omegaGrid, min = a, max = b)
  
  # calculate functions in the tangent space evaluated at the partition of [a, b]
  # the uniform distribution on [a, b] is used as the reference measure
  LFs <- matrix(0, nrow = N, ncol=2^(nd))
  for(i in 1:N){
    LFs[i, ] <- quantile(x = ObsList[[i]], pVals) - omegaGrid
  }
  
  # get the constraint matrix and vector
  A <- ConstraintMatrixGPCA(nd)
  b <- ConstraintVectorGPCA(a, b, nd)
  
  # perform a convex PCA and apply the k-means clustering to scores (initial clustering) 
  resultCPCAinitial <- ConvexPCA(M=initialM, X = LFs, A = A, b = b)
  initialClustering <- kmeans(resultCPCAinitial$scores, centers = K, algorithm = "MacQueen", nstart = nStart)
  clustConf0 <- as.factor(initialClustering$cluster) 
  indClustIds <- lapply(levels(clustConf0), function(u) which(clustConf0 == u))
  
  # perform CPCA for each cluster
  listOfCPCAobjs <- list()
  for(c in 1:K){
    listOfCPCAobjs[[c]] <- ConvexPCA(M = iterM, X = LFs[indClustIds[[c]], ], A, b)
  }
  
  # perform the reclassification
  clustConf <- list()
  iterNum <- 0
  for(j in 1:iterMax){
    
    # update the cluster memberships 
    iseCosts <- GetDistances(LFs, listOfCPCAobjs, A, b)
    clustConf[[j]] <- as.factor(apply(iseCosts, 1, which.min))
    indClustIds <- lapply(levels(clustConf[[j]]), function(u) which(clustConf[[j]] == u) )
    
    # update the cluster structures
    listOfCPCAobjs <- list()
    for(c in 1:K){
      listOfCPCAobjs[[c]] <- ConvexPCA(M = iterM, X = LFs[indClustIds[[c]], ], A, b)
    }
    
    iterNum <- iterNum + 1
    
    # check if algorithm converged
    curvesThatChanged <- ifelse(j > 1, sum(!( as.numeric(clustConf[[j]]) == as.numeric(clustConf[[j-1]] ))),
                                sum(!( as.numeric(clustConf[[j]]) == as.numeric(clustConf0))))
    if(curvesThatChanged == 0){
      break
    }
  }
  
  # store the result
  resultClustering <- list()
  resultClustering$cluster <- clustConf[[iterNum]]
  resultClustering$cpcaList <- listOfCPCAobjs
  resultClustering$iterToConv <- iterNum-1
  resultClustering$prevConf <- clustConf
  resultClustering$clustConf0 <- clustConf0
  return(resultClustering)
}















