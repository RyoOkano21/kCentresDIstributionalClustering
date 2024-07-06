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


kCDC <- function(ObsList, a, b, nd, K, initialM, iterM, nStart, iterMax){
  # Performs the proposed k centres distributional clustering
  
  # input
  # ObsList: A list of vectors containing the observations from each distribution
  # a: A scalar defining the left end point of the interval supporting distributions
  # b: A scalar defining the right end point of the interval supporting distributions
  # nd: A scalar defining the order of the dyadic partition of the interval supporting distributions
  # K: A scalar defining the number of clusters
  # initialM: A scalar defining the number of principal components used in the initial clustering 
  # iterM: A scalar defining the number of principal components used in the reclassification
  # nStart: A scalar defining the number of random sets chosen in the k-means clustering in the initial clustering
  # iterMax: A scalar defining the maximum number of iterations allowed
  
  # output
  # cluster: A vector of levels 1:k, indicating the cluster to which each distribution is allocated
  # cpcaList: A list with the CPCA structures (means, convex principal directions and so on) for each cluster
  # iterToConv: A number indicating how many iterations where required until convergence
  # prevConf: A list of vectors containing the results of the previous clusterings
  # clustConf0: A vector of levels 1:k, indicating the result of the initial clustering
  
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


DistributionCVE <- function(maxComp, ObsList, a, b, nd){
  # computes the cumulative variance explained 
  # maxComp is the maximum number of principal components requested
  # ObsList is a list of vectors containing the observations from each distribution
  
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