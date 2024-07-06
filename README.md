# _k_-Centres Distributional Clustering 
This repository provides R code implementing the _k_-centres distributioinal clustering, as proposed by the following paper. 

Okano, R. and Imaizumi, M. (2024). Wasserstein _k_-Centres Clustering for Distributional Data. 

The repository includes the following files.

- `HelperFunctions.cpp`: Helper functions for soloving the convex PCA problem with an algorithm proposed by [Campbell & Wong (2022)]
- `CPCAFunctions.R`: Functions for implenetations of convex PCA 
- `GPCAFunctions.R`: Functions for implenetations of geodesic PCA in the Wasserstein space
- `kCDC.R`: Functions for implementaitons of the proposed method
- `ClusteringFunctions.R`: Functions for implementations of some existing clustering methods
- `Simulation.R`: Example of applying the proposed and some existing clustering methods to a simulated data
- `RealDataAnalysis.R`: Example of applying the proposed and some existing clustering methods to a real data



