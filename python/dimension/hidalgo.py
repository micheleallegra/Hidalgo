import numpy as np
import random
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
import sys
#sys.path.append("./lib/python")
import _gibbs



class  hidalgo():

	def __init__(self, metric = 'euclidean', K = 2, zeta=0.8, q=3, Niter=10000, Nreplicas=1,burn_in=0.9,a=np.ones(2),b=np.ones(2),c=np.ones(2),f=np.ones(2)):

		a=np.ones(K)
		a=np.ones(K)
		a=np.ones(K)
		a=np.ones(K)

		self.metric = metric
		self.K = K
		self.zeta = zeta
		self.q = q
		self.Niter=Niter
		self.burn_in=burn_in
		self.Nreplicas=Nreplicas
		self.a=a
		self.b=b
		self.c=c
		self.f=f

	def _fit(self,X):

		q=self.q
		K=self.K
		Niter=self.Niter
		zeta=self.zeta
		Nreplicas=self.Nreplicas
		a=self.a
		b=self.b
		c=self.c
		f=self.f
		burn_in=self.burn_in

		assert isinstance(X,np.ndarray), "X should be a numpy array"
		assert len(np.shape(X))==2, "X should be a two-dimensional numpy array"

		N,d = np.shape(X)

		if self.metric!='predefined':
			nbrs = NearestNeighbors(n_neighbors=q+1, algorithm='ball_tree',metric=self.metric).fit(X)
			distances, indicesIn = nbrs.kneighbors(X)

			nbrmat=np.zeros((N,N))

			for n in range(q):
				nbrmat[indicesIn[:,0],indicesIn[:,n+1]]=1

			nbrcount=np.sum(nbrmat,axis=0)
			indicesOut=np.where(nbrmat.T)[1]
			indicesTrack=np.cumsum(nbrcount)
			indicesTrack=np.append(0,indicesTrack[:-1])
			
		else:
			distances = np.sort(X)[:,:q+1]
			indicesIn = np.argsort(X)[:,:q+1]

			nbrmat=np.zeros((N,N))

			for n in range(q):
				nbrmat[indicesIn[:,0],indicesIn[:,n+1]]=1

			nbrcount=np.sum(nbrmat,axis=0)
			indicesOut=np.where(nbrmat.T)[1]
			indicesTrack=np.cumsum(nbrcount)
			indicesTrack=np.append(0,indicesTrack[:-1])

		mu = np.divide(distances[:,2],distances[:,1])

		fixed_Z=0;
		use_Potts=1;
		estimate_zeta=0;
		sampling_rate=10;
		Nsamp=np.floor((Niter-np.ceil(burn_in*Niter))/sampling_rate).astype(int)
		self.Nsamp=Nsamp
		Npar=N+2*K+2+1*(estimate_zeta);

		sampling=2*np.ones(Nsamp*Npar);
		bestsampling=np.zeros((Nsamp,Npar));

		indicesIn=indicesIn[:,1:]
		indicesIn=np.reshape(indicesIn,(N*q,))

		maxlik=-1.E10

		for r in range(Nreplicas):
			_gibbs.GibbsSampling(Niter,K,fixed_Z,use_Potts,estimate_zeta,q,zeta,sampling_rate,burn_in,r,mu,indicesIn.astype(float),indicesOut.astype(float),nbrcount,indicesTrack,a,b,c,f,sampling)
			sampling=np.reshape(sampling,(Nsamp,Npar))	
			lik=np.mean(sampling[:,-1],axis=0)
			if(lik>maxlik): 
				bestsampling=sampling
			sampling=np.reshape(sampling,(Nsamp*Npar,))	

		return bestsampling



	def fit(self, X):
		N=np.shape(X)[0]
		sampling=self._fit(X)
		K=self.K

		self.d_=np.mean(sampling[:,:K],axis=0)
		self.derr_=np.std(sampling[:,:K],axis=0)
		self.p_=np.mean(sampling[:,K:2*K],axis=0)
		self.perr_=np.std(sampling[:,K:2*K],axis=0)
		self.lik_=np.mean(sampling[:,-1],axis=0)
		self.likerr_=np.std(sampling[:,-1],axis=0)
		
		Pi=np.zeros((K,N))

		for k in range(K):
			Pi[k,:]=np.sum(sampling[:,2*K:2*K+N]==k,axis=0)

		self.Pi=Pi/self.Nsamp
		Z=np.argmax(Pi,axis=0);
		pZ=np.max(Pi,axis=0);
		Z=Z+1
		Z[np.where(pZ<0.8)]=0
		self.Z=Z

		




'''
		if self.method=='MaxLikelihood':
			DimEstimate = float(N)/np.sum(np.log(mu))
		else:
			F = np.arange(1,N+1)/float(N)
			Neff = np.floor((1-self.discard)*N)
	
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
		Neff = np.floor((1-self.discard)*N)
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

'''

