// g++ -g -std=c++0x -o /scratch/efacco/ID_CALCULATION/NU_METHOD/BAYES/bayes_two_nn/Potts.x /scratch/efacco/ID_CALCULATION/NU_METHOD/BAYES/bayes_two_nn/Potts.cc

// Usage: ./Potts.x -input <filename> {-coord|-dist} [-discard <fraction>]

// Potts.cc: A K fissato implementiamo il metodo in 


#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include <functional>
#include <algorithm> 
#include <numeric>
#include <assert.h>
#include <limits>
#include "gibbs.cpp"
//#include <malloc.h>
#include <ctime>




using namespace std;


struct value_index{
        double value;
        int index;
};


bool customer_sorter_lower(value_index const& lhs, value_index const& rhs)
{
  return lhs.value < rhs.value;
}


bool customer_sorter_greater(value_index const& lhs, value_index const& rhs)
{
  return lhs.value > rhs.value;
}


void Gibbs(int K, vector <double> NU, double a[], double b[], double c[], vector <double> &sampling, int Niter);


double DistValue(int i, int j, int N, double *Dist, bool TriangSup)
{
	int ll, mm, kk;
	ll=max(i,j);
	mm=min(i, j);
	if(TriangSup) 
		kk=mm*(2*N-mm-3)/2.+ll;
	else
		kk=(ll*ll-ll)/2.+mm+1;
	return Dist[kk-1];
}


void find_nearest_dist(vector<double>& D, vector<int>& I, double *Dist, const int N, const int q, const bool TriangSup)
{
	for(int i=0; i<N; i++)
	{

		value_index dist_temp[N];  

		for(int j=0; j<N; j++)
		{
			double dist = DistValue(i,j,N,Dist,TriangSup); 
			if(j==i) dist=numeric_limits<double>::max();
			//cout << dist << " ";
			dist_temp[j].value=dist;
			dist_temp[j].index=j;		
		}

		partial_sort(dist_temp,dist_temp+q,dist_temp+N,customer_sorter_lower);

		for(int j=0; j< q; j++) 
		{
	    	D.push_back(dist_temp[j].value);
	    	I.push_back(dist_temp[j].index);
		}
	}
}



void find_nearest_coo(vector<double>& D,vector<int>& I, double *X, const int N, const int q, int ncoords, bool periodicB)
{

	double maxdist=numeric_limits<double>::max();;
	double L[ncoords];
	for(int cc=0; cc<ncoords; cc++)
	{
	    	L[cc]=2*M_PI;  
	}


	int i,j,cc;	
	double dist,Xtemp;

	const int CHUNK=1000;

	D.resize(q*N);
	I.resize(q*N);

	#pragma omp parallel default(none), shared(X,L,D,I,ncoords,maxdist,periodicB,cout), private(i,j,cc,dist,Xtemp)
	{ 
	#pragma omp for schedule (dynamic,CHUNK) 

	for(i=0; i<N; i++)
	{
		value_index dist_temp[N];


		for(j=0; j<N; j++)
		{
			if(j==i) 
			{ 
				dist=maxdist;	
				dist_temp[j].value=dist;				
				dist_temp[j].index=j;				
				continue;
			}

			dist=0.;
			for(cc=0; cc<ncoords; cc++)
			{
				Xtemp=0;
				Xtemp=X[i*ncoords+cc]-X[j*ncoords+cc];

				if(periodicB)
				{
					if(abs(Xtemp)>L[cc]*0.5) 
						if(X[i*ncoords+cc]>X[j*ncoords+cc]) Xtemp=L[cc]-Xtemp;
						else Xtemp=L[cc]+Xtemp;
				}
				dist+=Xtemp*Xtemp;
			}
			dist=sqrt(dist); 

			dist_temp[j].value=dist;				
			dist_temp[j].index=j;			
		}
	
		partial_sort(dist_temp,dist_temp+q,dist_temp+N,customer_sorter_lower);
		
		for(j=0; j< q; j++) 
		{
	    		D.at(q*i+j)=dist_temp[j].value;
	    		I.at(q*i+j)=dist_temp[j].index;
		}
	}

	} // end omp


}

//########################################################################################################################################

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {


	double *ctype=mxGetPr(prhs[0]);
	//int ctype1=(int)ctype[0];
        bool coordinates=(bool)ctype[0];

	
	double *X;
	int N;
	int ncoords;
	double *Dist;

	if(coordinates)  {
		mexPrintf("cc ");
		X=mxGetPr(prhs[1]);
		N=mxGetN(prhs[1]);
		ncoords=mxGetM(prhs[1]);
	}
	else {
		Dist=mxGetPr(prhs[1]);
		N=(int)((1.+sqrt(1.+8.*mxGetM(prhs[1])))*0.5);
	}
	
	mexPrintf("%d %d", N, ncoords);

	double *Kpar=mxGetPr(prhs[2]);
	int K=(int)Kpar[0];

	double *qpar=mxGetPr(prhs[3]);
	int q=(int)qpar[0];

	double *zetapar=mxGetPr(prhs[4]);
	double zeta=zetapar[0];

	double *Niterpar=mxGetPr(prhs[5]);
	int Niter=(int)Niterpar[0];

	double *Nreplicapar=mxGetPr(prhs[6]);
	int Nreplica=(int)Nreplicapar[0];



	plhs[0]= mxCreateDoubleMatrix(K,Nreplica,mxREAL);
	double *expected_d_sampling=mxGetPr(plhs[0]);

	plhs[1]= mxCreateDoubleMatrix(K,Nreplica,mxREAL);
	double *expected_p_sampling=mxGetPr(plhs[1]);

	plhs[2]= mxCreateDoubleMatrix(N*K,Nreplica,mxREAL);
	double *Z_sampling=mxGetPr(plhs[2]);

	plhs[3]= mxCreateDoubleMatrix(2,Nreplica,mxREAL);
	double *expected_L_sampling=mxGetPr(plhs[3]);



	double xx, yy;
	double dd;
	double rd;
	vector<double> Y;
        vector<double> D;
        vector<int> I;

	//int replica;

	bool periodicB=0;
	bool TriangSup=1;

	if(coordinates==true) find_nearest_coo( D,I, X, N, q, ncoords, periodicB);
	else  find_nearest_dist(D, I, Dist, N, q, TriangSup);

	
	mexPrintf("  Distance Matrix calculation done. \n");

	vector<int> Iout;

	for(int i=0;i <N; i++) 
	{
		for(int j=0;j <q; j++) 
		{
			Iout.push_back(-1);
		}
	}
	
	int count_neighbors[N];
	for(int i=0;i <N; i++) count_neighbors[i]=0;

	for(int i=0;i <N; i++)
        {
                for(int j=0;j <q; j++)
                {
                        int ind = I.at(q*i+j);
                        count_neighbors[ind]++;
                }
        }
 
   	int track_neighbors[N];
	int count=0;

        for(int i=0;i <N; i++) {
		track_neighbors[i]=count;
		count+=count_neighbors[i];
	}

        for(int i=0;i <N; i++) count_neighbors[i]=0;

	for(int i=0;i <N; i++) 
	{
		for(int j=0;j <q; j++) 
		{
			int ind = I.at(q*i+j); 
			Iout.at(track_neighbors[ind]+count_neighbors[ind])=i;
			count_neighbors[ind]++;
		}
	}


	// COMPUTING MU  
    
	vector<double> MU;

	double num, den;
	double mu;

	for( int i=0; i<N; i++)
	{


		num=D[q*i+1];
		den=D[q*i];

		mu=num/den;

		if(den==0) 	{ mexErrMsgTxt("there are identical points!! %");   }

		if(mu>=1.)	MU.push_back(mu);
	}



	// LOOP OVER REPLICAS *****************************************************************************************************

	for(int replica=0; replica < Nreplica; replica++) 
	{ 

		 
		// MONTE CARLO SAMPLING

		double b[K];
		double a[K];
		double c[K];
		double f[2];


		for(int k=0; k<K; k++) 
		{
			a[k]=1;
			b[k]=1;
			c[k]=1;
		}


		f[0]=1; 
		f[1]=1; 

		vector <double> sampling;

		bool fixed_Z=false;
		bool use_Potts=true;
		bool estimate_zeta=false;
		int sampling_rate=10;
		

		Gibbs_Potts(Niter, K, fixed_Z, use_Potts, estimate_zeta, q, MU, I, Iout, count_neighbors, track_neighbors, a, b, c, f, zeta, sampling, sampling_rate,replica);


		// computing expected valued

		


		for(int k=0; k<K; k++) expected_d_sampling[K*replica+k]=0;
		for(int k=0; k<K; k++) expected_p_sampling[K*replica+k]=0;
		for(int k=0; k<2; k++) expected_L_sampling[2*replica+k]=0;



		for(int i=0; i<N; i++) 
		{
			for(int k=0; k<K; k++) Z_sampling[(K*N)*replica+K*i+k]=0;
		}



		int L=sampling.size();


      		int Num_par;
                if(fixed_Z==false) Num_par=N+2*K+2;
                else Num_par=N+2*K+2;
		if(use_Potts==false) Num_par=N+2*K+2;
	

		int Niter_sampled = L/Num_par;


	
		int kept=0;
		for(int it=0; it<Niter_sampled; it++)
		{


		    int c1=0;
		    for(int k=0; k<K; k++)
			{
			if(it >=1) kept++;
			if(it >=1) expected_d_sampling[K*replica+k]+=sampling.at(Num_par*it+c1);
			c1++;
			}
			for(int k=0; k<K; k++)
			{
				if(it>=1) expected_p_sampling[K*replica+k]+=sampling.at(Num_par*it+c1);
			c1++;
			}

			for(int i=0; i<N; i++)
			{
				if(it>=1)
				{
					int Z0=sampling.at(Num_par*it+c1);
					Z_sampling[(N*K)*replica+K*i+Z0]++;
				}
			c1++;
			}

			for(int k=0; k<2; k++)
                        {
                                if(it>=1) expected_L_sampling[2*replica+k]+=sampling.at(Num_par*it+c1);
			c1++;
			}

		}


		for(int k=0; k<K; k++) expected_d_sampling[K*replica+k]=expected_d_sampling[K*replica+k]/(Niter_sampled-1);
		for(int k=0; k<K; k++) expected_p_sampling[K*replica+k]=expected_p_sampling[K*replica+k]/(Niter_sampled-1);
		for(int k=0; k<2; k++) expected_L_sampling[2*replica+k]=expected_L_sampling[2*replica+k]/(Niter_sampled-1);





		for(int i=0; i<N; i++) 
		{
			int mode;
			double val_mode=0;
			for(int k=0; k<K; k++) 
			{
				Z_sampling[(N*K)*replica+K*i+k]=Z_sampling[(N*K)*replica+K*i+k]/(Niter_sampled-1);
				if(Z_sampling[(N*K)*replica+K*i+k] > val_mode) 
				{
					val_mode=Z_sampling[K*i+k];
					mode=k;
				}

			}
		}

	}


}



