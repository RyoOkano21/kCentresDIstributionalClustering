# _k_-Centres Distributional Clustering 
This repository provides R code implementing the _k_-centres distributioinal clustering, as proposed by [2]. 

The repository includes the following files.
- `HelperFunctions.cpp`: Helper functions for implementations of an algorithm proposed by [1] for solving convex PCA
- `CPCAFunctions.R`: Functions for implenetations of convex PCA 
- `GPCAFunctions.R`: Functions for implenetations of geodesic PCA in the Wasserstein space
- `kCDC.R`: Functions for implementaitons of the proposed method
- `ClusteringFunctions.R`: Functions for implementations of some existing clustering methods
- `Simulation.R`: Example of applying the proposed and some existing clustering methods to a simulated data
- `RealDataAnalysis.R`: Example of applying the proposed and some existing clustering methods to a real data
- `PopulationData.RData`: Population dataset used in `RealDataAnalysis.R`

# References
[1] Campbell, S. and T.-K. L. Wang (2022). Effifient Convex PCA with application to Wasserstein geodesic PCA and ranked data. _arXiv preprint arXiv:2211.02990_.

[2] Okano, R. and M. Imaizumi (2024). Wasserstein _k_-Centres Clustering for Distributional Data. _arXiv preprint arXiv:????_.




