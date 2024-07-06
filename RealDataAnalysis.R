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
load("PopulationData.RData")


### Process data ###

# Get observations from each distribution 
ObsListMale <- list()
for(i in 1:42){
  ageVec <- c()
  for(j in 1:101){
    ageVec <- c(ageVec, rep(j-1, PopulationMale[i, j+1]))
  }
  ObsListMale[[i]] <- ageVec
}

ObsListFemale <- list()
for(i in 1:42){
  ageVec <- c()
  for(j in 1:101){
    ageVec <- c(ageVec, rep(j-1, PopulationFemale[i, j+1]))
  }
  ObsListFemale[[i]] <- ageVec
}

# Plot distributions
for(i in 1:42){
  plot(density(ObsListMale[[i]], adjust = 1, from = 0, to = 100), xlim=c(0, 100), 
       ylim=c(0, 0.025), main="", xlab="", ylab="", xaxt="n", yaxt="n", bty="n");par(new=T)
}
for(i in 1:41){
  plot(density(ObsListFemale[[i]], adjust = 1, from = 0, to = 100), xlim=c(0, 100), 
       ylim=c(0, 0.025), main="", xlab="", ylab="", xaxt="n", yaxt="n", bty="n", col="red", lty=2);par(new=T)
}
plot(density(ObsListFemale[[42]], adjust = 1, from = 0, to = 100), xlim=c(0, 100), 
     ylim=c(0, 0.025), main="", xlab="age", ylab="density", col="red", lty=2)
legend("topright", legend = c("men", "women"), col = c("black", "red"), 
       lty = c(1, 2), pt.cex = 1.5, cex=1.5)

# Combine the data set
ObsList <- c(ObsListMale, ObsListFemale)


### Perform clustering ###

# Choose the number of principal components based on the cumulative variation explained 
tau <- 0.85
CVE <- DistributionCVE(maxComp = 10, ObsList = ObsList, a = 0, b = 100, nd = 4)
M <-  which(CVE >= tau)[1] 

# Perform clustering using the proposed method
resultKCDC <- kCDC(ObsList, a = 0, b = 100, nd = 4, K = 2, initialM = M, iterM = M, 
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
labelTrue <- c(rep(1, 42), rep(2, 42))

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


### Visualize the estimated structure of each cluster  ###

# get the mean functions and convex principal directions in the tangent space 
# evaluated at the partition of [0, 100]
meanClust1 <- resultKCDC$cpcaList[[1]]$mean
meanClust1 <- c(0, meanClust1)
meanCLust2 <- resultKCDC$cpcaList[[2]]$mean
meanClust2 <- c(0, meanCLust2)
CPdirecClust1 <- resultKCDC$cpcaList[[1]]$directions
CPdirecClust1 <- c(0, CPdirecClust1)
CPdirecClust2 <- resultKCDC$cpcaList[[2]]$directions
CPdirecClust2 <- c(0, CPdirecClust2)

# define functions by adding the identity function to the mean functions 
OmegaGrid <- DyadicPartition(a = 0, b = 100, n = 4)
meanIdClust1 <- function(x){
  interp1(x = OmegaGrid, y = meanClust1+OmegaGrid, xi = x)
}

meanIdClust2 <- function(x){
  interp1(x = OmegaGrid, y = meanClust2+OmegaGrid, xi = x)
}

# generate samples from the Fr\'{e}chet mean of each cluster
set.seed(123)
sampleUnif <- runif(n=10^5, min = 0, max = 100)
sampleMeanClust1 <- sapply(sampleUnif, meanIdClust1)
sampleMeanClust2 <- sapply(sampleUnif, meanIdClust2)

# plot the densities of the Fr\'{e}chet mean of each cluster
plot(density(sampleMeanClust1, adjust = 2, from = 0, to = 100), 
     ylim=c(0, 0.025), main="", xlab="", ylab="", cex.axis=1, lwd=2);par(new=T)
plot(density(sampleMeanClust2, adjust = 2, from = 0, to = 100), ylim=c(0, 0.025), 
     col="red", main="", xlab="age", ylab="density", cex.axis=1, cex.lab=1, lwd=2, lty=2)
legend("topright", legend = c("Cluster 1", "Cluster 2"), col = c(1, 2), lty=c(1, 2), 
       lwd=c(2,2), cex = 1.25)

# define functions needed to plot the (fisrt) mode of variation of each cluster
sdClust1 <- sqrt(resultKCDC$cpcaList[[1]]$variances)
sdClust2 <- sqrt(resultKCDC$cpcaList[[2]]$variances)
CPdirecIdClust1Plus <- function(x){
  interp1(x = OmegaGrid, y = meanClust1+sdClust1*CPdirecClust1+OmegaGrid, xi = x)
}
CPdirecIdClust1Minus <- function(x){
  interp1(x = OmegaGrid, y = meanClust1-sdClust1*CPdirecClust1+OmegaGrid, xi = x)
}
CPdirecIdClust2Plus <- function(x){
  interp1(x = OmegaGrid, y = meanClust2+sdClust2*CPdirecClust2+OmegaGrid, xi = x)
}
CPdirecIdClust2Minus <- function(x){
  interp1(x = OmegaGrid, y = meanClust2-sdClust2*CPdirecClust2+OmegaGrid, xi = x)
}

# generate samples to plot the mode of variation of each cluster
sampleModeClust1Plus <- sapply(sampleUnif, CPdirecIdClust1Plus)
sampleModeClust1Minus <- sapply(sampleUnif, CPdirecIdClust1Minus)
sampleModeClust2Plus <- sapply(sampleUnif, CPdirecIdClust2Plus)
sampleModeClust2Minus <- sapply(sampleUnif, CPdirecIdClust2Minus)

# plot the mode of variation of each cluster
plot(density(sampleMeanClust1, adjust = 2, from = 0, to = 100), ylim=c(0, 0.025), main="", 
     lwd=2, xlab="", ylab="", xaxt="n", yaxt="n");par(new=T)
plot(density(sampleModeClust1Plus, adjust = 2, from = 0, to = 100), ylim=c(0, 0.025), main="", 
     lwd=1, xlab="", ylab="", xaxt="n", yaxt="n");par(new=T)
plot(density(sampleModeClust1Minus, adjust = 2, from = 0, to = 100), ylim=c(0, 0.025), 
     lwd=1, main="Cluster 1", xlab="age", ylab="density")
x_seg <- density(x = sampleMeanClust1, adjust = 2, from = 0, to = 100)$x
y1_seg <- density(x = sampleMeanClust1, adjust = 2, from = 0, to = 100)$y
y2_seg <- density(x = sampleModeClust1Plus, adjust = 2, from = 0, to = 100)$y
y3_seg <- density(x = sampleModeClust1Minus, adjust = 2, from = 0, to = 100)$y
segments(x_seg, y1_seg, x_seg, y2_seg, lwd = 0.3, lty = 1)
segments(x_seg, y1_seg, x_seg, y3_seg, lwd=0.3, lty = 1) 

# plot the mode of variation of each cluster
plot(density(sampleMeanClust2, adjust = 2, from = 0, to = 100), ylim=c(0, 0.025), main="", 
     lwd=2, xlab="", ylab="", xaxt="n", yaxt="n");par(new=T)
plot(density(sampleModeClust2Plus, adjust = 2, from = 0, to = 100), ylim=c(0, 0.025), main="", 
     lwd=1, xlab="", ylab="", xaxt="n", yaxt="n");par(new=T)
plot(density(sampleModeClust2Minus, adjust = 2, from = 0, to = 100), ylim=c(0, 0.025), 
     lwd=1, main="Cluster 2", xlab="age", ylab="density")
x_seg <- density(x = sampleMeanClust2, adjust = 2, from = 0, to = 100)$x
y1_seg <- density(x = sampleMeanClust2, adjust = 2, from = 0, to = 100)$y
y2_seg <- density(x = sampleModeClust2Plus, adjust = 2, from = 0, to = 100)$y
y3_seg <- density(x = sampleModeClust2Minus, adjust = 2, from = 0, to = 100)$y
segments(x_seg, y1_seg, x_seg, y2_seg, lwd = 0.3, lty = 1)
segments(x_seg, y1_seg, x_seg, y3_seg, lwd=0.3, lty = 1) 