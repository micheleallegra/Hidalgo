import numpy as np
from dimension import TwoNN
from dimension import hidalgo
from sklearn.metrics.pairwise import euclidean_distances

N=500

d1=1
d2=3

X=np.zeros((2*N,10))

for j in range(d1):
	X[:N,j]= np.random.normal(0,1,N)

for j in range(d2):
	X[N:,j]= np.random.normal(2,1,N)


K=2;

########## run TwoNN ###################################################################

model = TwoNN()

m = TwoNN(metric = 'euclidean', method = 'LinearFit',discard=0.1,block_analysis=False)

model.fit(X[:N,:])

print model.DimEstimate_

model.ShowLinearFit()

model.ShowBlockAnalysis()


############## run Hidalgo ##################################################################


model=hidalgo(K=K)

#model=hidalgo(K=2,Niter=2000,zeta=0.65,q=5,Nreplicas=10,burn_in=0.8)


model.fit(X)

print model.d_,model.derr_
print model.p_,model.perr_
print model.lik_, model.likerr_
print model.Pi
print model.Z


D = euclidean_distances(X)

model=hidalgo(metric = 'predefined',K=K)

model.fit(D)

print model.d_,model.derr_
print model.p_,model.perr_
print model.lik_, model.likerr_
print model.Pi
print model.Z

