# Reset Environment
rm(list=ls())

# Import source files
library(Rcpp)
library(truncnorm)
library(mclust)
sourceCpp("HelperFunctions.cpp")
source("CPCAFunctions.R")
source("GPCAFunctions.R")
source("kCDC.R")
source("ClusteringFunctions.R")

# Set seed
set.seed(123)



### Create data set ###

# Create data set of cluster 1
m <- function(x){
  qtruncnorm(x, a=0, b=1, mean = 0.75, sd = 0.3) - x
} 
phi1  <-  function(x){ 
  sqrt(2) * sin(2*pi*x)
}
phi2 <- function(x){
  sqrt(2) * sin(8*pi*x)
}
c1 <- 0.04
c2 <- 0.001
sampleClust1 <- GenerateDistributions(numSample = 50, m, phi1, phi2, 
                                      c1, c2, a = 0, b = 1, numObs = 2000)
# Create data set of cluster 2
m <- function(x){
  qtruncnorm(x, a=0, b=1, mean = 0.75, sd = 0.25) - x
} 
phi1  <-  function(x){ 
  sqrt(2) * sin(4*pi*x)
}
phi2 <- function(x){
  sqrt(2) * sin(6*pi*x)
}
c1 <- 0.02
c2 <- 0.0133
sampleClust2 <- GenerateDistributions(numSample = 50, m, phi1, phi2, 
                                      c1, c2, a = 0, b = 1, numObs = 2000)

# Plot distributions
for(i in 1:49){
  plot(density(sampleClust1[[i]], adjust = 2, from = 0, to = 1), xlim=c(0, 1), ylim=c(0, 3), main = "", 
       xlab = "", ylab = "", xaxt="n", yaxt="n", bty="n");par(new=T)
}
plot(density(sampleClust1[[50]], adjust = 2, from = 0, to = 1), xlim=c(0, 1), ylim=c(0, 3), 
     main = "Cluster 1", xlab = "", ylab = "", cex = F)
for(i in 1:49){
  plot(density(sampleClust2[[i]], adjust = 2, from = 0, to = 1), xlim=c(0, 1), ylim=c(0, 3), main = "", 
       xlab = "", ylab = "", xaxt="n", yaxt="n", bty="n");par(new=T)
}
plot(density(sampleClust2[[50]], adjust = 2, from = 0, to = 1), xlim=c(0, 1), ylim=c(0, 3), 
     main = "Cluster 2", xlab = "", ylab = "", cex = F)

# Combine the data set
ObsList <- c(sampleClust1, sampleClust2)


### Perform clustering ###

# Choose the dimension of models based on the cumulative variation explained
tau <- 0.9
CVE <- DistributionCVE(maxComp = 10, ObsList = ObsList, a = 0, b = 1, nd = 4)
M <-  which(CVE >= tau)[1] 

# Perform clustering using the proposed method
resultKCDC <- kCDC(ObsList, a = 0, b = 1, nd = 4, K = 2, initialM = M, iterM = M, 
                   nStart = 20, iterMax = 20)

# Perform clustering using other methods
resultWkmeans <- Wkmeans(ObsList, a = 0, b = 1, nd = 4, K = 2, nStart = 20, iterMax = 20) # Wassetstein k-means
resultTrim1 <- WkmeansTrim(ObsList, a = 0, b = 1, nd = 4, 
                           K = 2, nStart = 20, iterMax = 20, delta = 0.01) # trimmed Wasserstein k-means with delta = 0.01
resultTrim2 <- WkmeansTrim(ObsList, a = 0, b = 1, nd = 4, 
                           K = 2, nStart = 20, iterMax = 20, delta = 0.05) # trimmed Wasserstein k-means with delta = 0.05
resultTrim3 <- WkmeansTrim(ObsList, a = 0, b = 1, nd = 4, 
                           K = 2, nStart = 20, iterMax = 20, delta = 0.1) # trimmed Wasserstein k-means with delta = 0.1


### Compare the results ###

# define the true labels
labelTrue <- c(rep(1, 50), rep(2, 50))

# calculate the correct classification rates
1 - classError(labelTrue, resultKCDC$clustConf0)$errorRate # the initial clustering
1 - classError(labelTrue, resultKCDC$cluster)$errorRate # the proposed method
1 - classError(labelTrue, resultWkmeans)$errorRate # Wassetstein k-means
1 - classError(labelTrue, resultTrim1)$errorRate # trimmed Wasserstein k-means with delta = 0.01
1 - classError(labelTrue, resultTrim2)$errorRate # trimmed Wasserstein k-means with delta = 0.05
1 - classError(labelTrue, resultTrim3)$errorRate # trimmed Wasserstein k-means with delta = 0.1

# calculate the adjusted Rand index
adjustedRandIndex(labelTrue, resultKCDC$clustConf0) # the initial clustering
adjustedRandIndex(labelTrue, resultKCDC$cluster) # the proposed method
adjustedRandIndex(labelTrue, resultWkmeans) # Wassetstein k-means
adjustedRandIndex(labelTrue, resultTrim1) # trimmed Wasserstein k-means with delta = 0.01
adjustedRandIndex(labelTrue, resultTrim2) # trimmed Wasserstein k-means with delta = 0.05
adjustedRandIndex(labelTrue, resultTrim3) # trimmed Wasserstein k-means with delta = 0.1