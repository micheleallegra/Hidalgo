# Hidalgo
find clusters of different intrinsic dimension
https://arxiv.org/abs/1902.10459


COMPILE SOURCE
* mex -v CXX='$CXX -fopenmp' CXXOPTIMFLAGS='-O3 -DNDEBUG' hidalgo.cpp

SYNTAX

 
* out=HeterogeneousID(X): apply Hidalgo to coordinate matrix X (rows=observations, cols=coordinates) 
* out=HeterogeneousID('Distance',D): apply Hidalgo to distance matrix D  
* out=HeterogeneousID(X,'K',K): apply Hidalgo to coordinate matrix X, number of clusters is set to K (by default an optimal K is searched) 
* out=HeterogeneousID(X,'q',Q): apply Hidalgo to coordinate matrix X, number of neighbors is set to q (by default q=3) 
* out=HeterogeneousID(X,'Zeta',z): apply Hidalgo to coordinate matrix X, neighborhood uniformity is set to z (by default z=0.8)
* out=HeterogeneousID(X,'Niter',N): apply Hidalgo to coordinate matrix X, numer of Gibbs iterations is set to N (by default N=10,000) 
* out=HeterogeneousID(X,'Nreplicas',R): apply Hidalgo to coordinate matrix X, repeating the clustering Nreplicas times and keeping the optimal result (by default Nreplicas=10)

OUTPUT
* out.d: intrinsic dimensions of the clusters
* out.p: sizes of the cluster
* out.Z: probability of assignment of each point to the clusters (rows=clusters,cols=points)
