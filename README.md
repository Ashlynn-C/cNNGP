# cNNGP

* The *CNNGP* folder started with the [spNNGP](https://cran.r-project.org/package=spNNGP) package and has been edited to include the code for the clustered NNGP.
* The *leaderCluster2* folder has code from the [leaderCluster](https://cran.r-project.org/package=leaderCluster) package with some very minor edits.
* The *GEDI* folder contains the data preparation, model fitting, and summarizing code for our analysis of the GEDI biomass data.
* The *experiments* folder contains the code for the simulated data generation, model fitting, and analysis. 
* The *testing.R* file contains example code to get started using the *clustNNGP* function. 

All of the original functions in the *spNNGP* package still work in the *CNNGP* package. Additionally, there are two new functions available:

* *clustering*
* *clustNNGP*

***Important***

Do not order the coords before passing to *clustNNGP*. The order must be specified as the *ord* parameter. If *ord* is not given, *clustNNGP* will run but will order the points by  

    ord <- order(coords[,1])##default order by x column

## *clustering*

This function takes the coordinates and returns the approximated nearest neighbor distances. It requires the following arguments:
* **coords**: Coordinates of the data. Can be ordered or ordering vector can be passed as argument. 
* **ord** (Optional): Coordinate ordering vector. 
* **m**: Number of nearest neighbors. 
* **indx** (Optional): If the nearest neighbors structure has been precomputed, it can be passed as an argument here.
* **radius**: Radius value for leader clustering algorithm.
* **clust.method**: Can be "clustering" or "exact". If exact, the original distances are used instead of clusters. Default is "clustering".
* **PC**: Boolean. Whether the clustering should be done on the principal components or original distances. Default is TRUE.
* **PCvar**: Percent of variance in original distances explained by the principal components. Default is 0.9.
* **max.iter**: Number of iterations for the leader clustering algorithm. Default is 1.
* **search.type**: Nearest neighbor searching algorithm. Default is same as the default in the clustNNGP calls.
* **omp.threads** (Optional): Number of threads to use to compute nearest neighbor index structure. Default is 1. 

This function returns

* **cluster_id**: The cluster assignment of each row of coords.
* **cluster_idx**: The row indices of the cluster centroids. Used to get centroids in terms of the original distances since clustering can be done on the principal components. 
* **cluster_centroids**: Centroids from the clustering algorithm. 
* **num_clusters**: Number of cluster centroids needed in the clustering. 
* **iter**: Number of clustering iterations plus 1. 

## *clustNNGP*

This version of the NNGP only works for a sequential model for the Gaussian family. New arguments are:

* **clust.method**: Can be "clustering" or "exact". If exact, the original distances are used instead of clusters. Default is "clustering".
* **radius**: Radius value for leader clustering algorithm.
* **PC**: Boolean. Whether the clustering should be done on the principal components or original distances. Default is TRUE.
* **PCvar**: Percent of variance in original distances explained by the principal components. Default is 0.9.
* **max.iter**: Number of iterations for the leader clustering algorithm. Default is 1.
* **clustering.info**: If output from *clustering* function is precomputed, it can be passed to the clustNNGP function here. 