library(stats)

Wkmeans <- function(ObsList, a, b, nd, K, nStart, iterMax){
  # Performs the Wasserstein k-means clustering
  
  # input 
  # ObsList: A list of vectors containing the observations from each distribution
  # a:  A scalar defining the left end point of the interval supporting distributions
  # b: A scalar defining the right end point of the interval supporting distributions
  # nd: A scalar defining the order of the dyadic partition of the interval supporting distributions
  # K: A scalar defining the number of clusters
  # nStart: A scalar defining the number of random sets chosen
  # iterMax: A scalar defining the maximum number of iterations allowed
  
  # output: A vector of levels 1:K, indicating the cluster to which each distribution is allocated
  
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
  # ObsList: A list of vectors containing the observations from each distribution
  # a:  A scalar defining the left end point of the interval supporting distributions
  # b: A scalar defining the right end point of the interval supporting distributions
  # nd: A scalar defining the order of the dyadic partition of the interval supporting distributions
  # K: A scalar defining the number of clusters
  # nStart: A scalar defining the number of random sets chosen
  # iterMax: A scalar defining the maximum number of iterations allowed
  # delta: A scalar defining the trimming level
  
  # output: A vector of levels 1:K, indicating the cluster to which each distribution is allocated
  
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


GenerateDistributions <- function(numSample, m, phi, lambda, a, b, numObs){
  # Given the structure of cluster, this function first generates distributions
  # belonging to the cluster, and then generates samples from the generated distributions
  
  # input
  # numSample: A scalar defining the number of distributions generated
  # m: A function defining the mean structure of cluster in the tangent space
  #    (the uniform distribution on the interval [a, b] is used as the reference measure)
  # phi: A function defining the convex principal directions of the cluster
  # lambda: A vector defining the ranges of the coefficients
  # a: A scalar defining the left end point of the interval supporting distributions
  # b: A scalar defining the right end point of the interval supporting distributions
  # numObs: A scalar defining the size of samples from each distribution
  
  # output: A list of vectors containing observations from each generated distribution
  
  J <- length(lambda)
  listSample <- list()
  for(i in 1:numSample){
    # generate the coefficients
    xi <- numeric(J)
    for(j in 1:J){
      xi[j] <- runif(n=1, min = -lambda[j], max = lambda[j])
    }
    
    # set the function in the tangent space
    g <- function(x){
      vals <- numeric(J)
      for(j in 1:J){
        vals[j] <- phi(x, j)
      }
      m(x) + xi　%*% vals
    }
    
    # generate the sample 
    sampleUnif <-  runif(n=numObs, min = a, max = b)
    listSample[[i]] <- sapply(sampleUnif, function(x){g(x)+x})
  }
  return(listSample)
} 