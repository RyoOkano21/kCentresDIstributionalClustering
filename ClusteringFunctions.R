library(stats)

Wkmeans <- function(ObsList, a, b, nd, K, nStart, iterMax){
  # Performs the Wasserstein k-means clustering
  
  # input 
  # ObsList: 
  # a: 
  # b: 
  # nd: 
  # K: 
  # nStart: 
  # iterMax
  
  N <- length(ObsList)
  
  # get the Dyadic partitions of the intervals [a, b] and [0, 1]
  omegaGrid <- DyadicPartition(a, b, nd)[-1]
  pVals <- punif(q = omegaGrid, min = a, max = b)
  
  # calculate the quantile functions evaluated at the partition of [0, 1]
  QFs <- matrix(0, nrow = N, ncol=2^(nd))
  for(i in 1:N){
    QFs[i, ] <- quantile(x = ObsList[[i]], pVals)
  }
  
  # perform the k-means clustering
  resultKmeans <- kmeans(x = QFs, centers = K, iter.max = iterMax, 
                         nstart = nStart, algorithm = "MacQueen")
  return(resultKmeans$cluster)
}


WkmeansTrim <- function(ObsList, a, b, nd, K, nStart, iterMax, delta){
  # Performs the k-means clustering with the trimmed Wasserstein distance
  
  # input 
  # ObsList: 
  # a: 
  # b: 
  # nd: 
  # K: 
  # nStart: 
  # iterMax: 
  # delta: 
  
  N <- length(ObsList)
  
  # get the Dyadic partitions of the intervals [a, b] and [0, 1]
  omegaGrid <- DyadicPartition(a, b, nd)[-1]
  pVals <- punif(q = omegaGrid, min = a, max = b)
  
  # calculate the quantile functions evaluated at the partition of [delta, 1-delta]
  condTrim <- (pVals >= delta) & (pVals <= (1-delta))
  numTrim <- sum(condTrim)
  QFsTrim <- matrix(0, nrow = N, ncol = numTrim)
  for(i in 1:N){
    QFsTrim[i, ] <- quantile(x=ObsList[[i]], pVals)[condTrim]
  }
  
  # perform the k-means clustering
  resultKmeansTrim <- kmeans(x = QFsTrim, centers = K, iter.max = iterMax, 
                             nstart = nStart, algorithm = "MacQueen")
  return(resultKmeansTrim$cluster)
}


GenerateDistributions <- function(numSample, m, phi1, phi2, c1, c2, a, b, numObs){
  # generates
  
  listSample <- list()
  for(i in 1:numSample){
    # generate the coefficients
    xi1 <- runif(n=1, min=-c1, max = c1)
    xi2 <- runif(n=1, min=-c2, max = c2)
    
    # set the function in the tangent space
    g <- function(x){
      m(x) + xi1*phi1(x) + xi2*phi2(x)
    }
    
    # generate the sample 
    sampleUnif <-  runif(n=numObs, min = a, max = b)
    listSample[[i]] <- sapply(sampleUnif, function(x){g(x)+x})
  }
  return(listSample)
} 
