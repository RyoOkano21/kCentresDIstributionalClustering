# _k_-Centres Distributional Clustering 
This repository provides R codes implementing the _k_-centres distributioinal clustering, as proposed by [2]. 

The repository includes the following files.
- `HelperFunctions.cpp`: Helper functions used in `CPCAFunctions.R` 
- `CPCAFunctions.R`: Functions for implenetations of convex PCA
- `GPCAFunctions.R`: Functions for implenetations of geodesic PCA in the Wasserstein space
- `kCDC.R`: Functions for implementaitons of the proposed method
- `ClusteringFunctions.R`: Functions for implementations of some existing clustering methods
- `Simulation.R`: Example of applying the proposed and some existing clustering methods to a simulated data
- `RealDataAnalysis.R`: Example of applying the proposed and some existing clustering methods to a real data
- `PopulationData.RData`: Population dataset used in `RealDataAnalysis.R`

For solving convex PCA problems, we use an algorithm developed by [1]. 
The dataset in `PopulationData.RData` is downloaded from STATcube – Statistical Database of Statistics Austria (https://www.statistik.at/en/databases/statcube-statistical-database).

Please contact `okano-ryo1134 (at) g.ecc.u-tokyo.ac.jp` if you have any comments.

# References
[1] Campbell, S. and T.-K. L. Wang (2024). Effifient Convex PCA with application to Wasserstein geodesic PCA and ranked data. _Journal of Computational and Graphical Statistics, 34_(2), 540–551.

[2] Okano, R. and M. Imaizumi (2025). Wasserstein _k_-Centres Clustering for Distributional Data. _Statistics and Computing, 35_.




