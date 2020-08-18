import numpy as np
import random
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt


class  TwoNN():

	def __init__(self, metric = 'euclidean', method = 'LinearFit',discard=0.2,block_analysis=True):
		self.metric = metric
		self.method = method
		self.block_analysis=block_analysis
		self.discard = discard

	def _fit(self,X):
		assert isinstance(X,np.ndarray), "X should be a numpy array"
		assert len(np.shape(X))==2, "X should be a two-dimensional numpy array"

		N,d = np.shape(X)

		if self.metric!='predefined':
			nbrs = NearestNeighbors(n_neighbors=3, algorithm='ball_tree',metric=self.metric).fit(X)
			distances, indices = nbrs.kneighbors(X)
		else:
			distances = np.sort(X)[1:3]
			indices = np.argsort(X)[1:3]

		mu = np.divide(distances[:,2],distances[:,1])
		mu = np.sort(mu)

		if self.method=='MaxLikelihood':
			DimEstimate = float(N)/np.sum(np.log(mu))
		else:
			F = np.arange(1,N+1)/float(N)
			Neff = np.floor((1-self.discard)*N).astype(int)
	
			par = np.polyfit(np.log(mu[:Neff]),-np.log(1-F[:Neff]),1)
			DimEstimate = par[0]

		return  DimEstimate,mu


	def fit(self, X):
		DimEstimate, mu= self._fit(X)
		self.DimEstimate_ = DimEstimate
		self.mu_ = mu
		N,d = np.shape(X)


		
	
		if self.block_analysis==True: 
			maxnblocks=20
			BlockSize =  np.zeros(maxnblocks)
			BlockDimMeanEstimate = np.zeros(maxnblocks)
			BlockDimStdEstimate = np.zeros(maxnblocks)

			idx = range(N)
			random.shuffle(idx)

			for nblocks in range(1,maxnblocks+1):
				BlockDimEst=[]
				blocksize=N/nblocks

				for b in range(nblocks):
					Y=X[idx[b*blocksize:(b+1)*blocksize],:]
					bde,mu=self._fit(Y) 					
					BlockDimEst=np.append(BlockDimEst,bde)

				BlockSize[nblocks-1]=blocksize
				BlockDimMeanEstimate[nblocks-1]=np.mean(BlockDimEst)
				BlockDimStdEstimate[nblocks-1]=np.std(BlockDimEst)	
				
			self.BlockEstimates_=np.column_stack((BlockSize,BlockDimMeanEstimate,BlockDimStdEstimate))

		return self


	def ShowLinearFit(self):
		mu=self.mu_
		DimEstimate=self.DimEstimate_
		N = np.shape(mu)[0]
		plt.figure()		
		F = np.arange(1,N+1)/float(N)
		Neff = np.floor((1-self.discard)*N).astype(int)
		plt.xlabel("log(mu)")
		plt.ylabel("log(1-F)")
		plt.title("d = %(DimEstimate)2.2f" % {"DimEstimate": DimEstimate})
		plt.scatter(np.log(mu[:Neff]),-np.log(1-F[:Neff]),color="red")
		plt.scatter(np.log(mu[Neff:-1]),-np.log(1-F[Neff:-1]),color="gray")
		plt.show()
		
	def ShowBlockAnalysis(self):
			BlockEstimates=self.BlockEstimates_
			plt.figure(2)		
			plt.xlabel("block size")
			plt.xscale('log')
			plt.ylabel("d estimate")
			plt.grid()
			plt.plot(BlockEstimates[:,0],np.ones(np.shape(BlockEstimates)[0])*BlockEstimates[0,1],color="blue")
			plt.errorbar(BlockEstimates[:,0],BlockEstimates[:,1],BlockEstimates[:,2],color="red")

			plt.show()		


