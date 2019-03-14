#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include <functional>
#include <algorithm>
#include <numeric>
#include <assert.h>
#include "mex.h"

using namespace std;


// binomial coefficient
double binom(int N, int q) {

 double s=1;
 
 if(q==0) return s;
 
 else{ 
 	for(int q1=0; q1<q; q1++) {
              s=s*(N-q1)/(q1+1);
   	}
	return s; 
 }

}


// partition function Z
double Zpart(int N, int N1, double zeta, int q){

        double s=0;

        for(int q1=0; q1<=q; q1++) {
	     s=s+binom(N1-1,q1)*binom(N-N1,q-q1)*pow(zeta,q1)*pow(1-zeta,q-q1);	

        }
        return s;
}


/*

K ---------------------------------------: number of manifolds
fixed_Z ---------------------------------: estimate parameters with fixed Z
use_Potts -------------------------------: use local interaction between Z
estimate_zeta ---------------------------: update z in the sampling
q ---------------------------------------: number of points for local Z interaction
MU --------------------------------------: r2/r1
Iin -------------------------------------: matrix of neighbors of point i
Iout ------------------------------------: matrix of points for which i is neighbor
Iout_count-------------------------------: matrix of points for which i is neighbor (index tracking)
Iout_track-------------------------------: matrix of points for which i is neighbor (index tracking)
a, b ------------------------------------: prior parameters of d 
c ---------------------------------------: prior parameters of p
f ---------------------------------------: parameters of zeta
zeta ------------------------------------: fixed zeta value
sampling --------------------------------: sampling record
sampling_rate ---------------------------: sampling record rate
*/


void Gibbs_Potts(int Niter, int K, bool fixed_Z, bool use_Potts, bool estimate_zeta, int q, vector <double> MU, vector<int> Iin, vector<int> Iout, int Iout_count[], int Iout_track[],  double a[], double b[], double c[], double f[], double zeta, vector <double>& sampling, int sampling_rate, int replica){


	srand(0);

	int N=MU.size();

	double d[K];
	double V[K];	
	int NN[K];
	int Z[N];
	double a1[K];
	double b1[K];
	double c1[K];
	double p[K];
	double pp;
	double f1[2]; 




/*************** INITIALIZE PARAMETERS *************************************************************************************/

	srand (10000*replica);

	for(int k=0; k<K; k++) { 
		V[k]=0;  
		NN[k]=0;
		d[k]=1;
		sampling.push_back(d[k]);
	}

	for(int k=0; k<K; k++) {
		p[k]=1./K;
		sampling.push_back(p[k]);
	}

	pp=(double(K-1))/K;

  //  if(use_Potts==true) sampling.push_back(zeta);

	
	for(int i=0; i<N; i++){
		if(fixed_Z==false) { 
			int z=rand()%K;
			Z[i]=z;
		}
		else Z[i]=0;
		//cout << Z[i] << endl;;

		V[Z[i]]=V[Z[i]]+log(MU[i]);
		NN[Z[i]]+=1;
		sampling.push_back(Z[i]);
	}

	for(int k=0; k<K; k++) {
		a1[k]=a[k]+NN[k];
		b1[k]=b[k]+V[k];
		c1[k]=c[k]+NN[k];
	}

	double N_in=0;
	for(int i=0; i< N; i++) {
		int k=Z[i];
		for(int j=0; j<q;j++) {
			int index=Iin.at(q*i+j);
			if(Z[index]==k) N_in=N_in++;
	    }			
	}	


	f1[0]=f[0]+N_in; 
	f1[1]=f[1]+N*q-N_in;

        sampling.push_back(0);
        sampling.push_back(0);


/******************* ITERATIONS **********************************************************************************************************/
	
	for(int it=0; it<Niter; it++){

		if(it%100==0) mexPrintf("it= %d \n",it);

		bool stop;


		/* SAMPLING d***********************************************************************************************************/

		for(int k=0; k<K; k++){
	        
			stop=false;


			while(stop==false){
	   
				double r1 = double(rand())/RAND_MAX*200;
				double r2 = double(rand())/RAND_MAX;
	      
				double rmax = (a1[k]-1)/b1[k];
				double frac;

				if(a1[k]-1>0) frac = exp(-b1[k]*(r1-rmax)-(a1[k]-1)*(log(rmax)-log(r1))); 
				else frac = exp(-b1[k]*r1);


				if(frac>r2){  
					stop=true;
					if(it%sampling_rate==0 &&  it>= Niter/10*9) sampling.push_back(r1);
					d[k]=r1; 
				}
			}
		}


		/* SAMPLING p **********************************************************************************************************/


		for(int k=0; k<K-1; k++) {

			stop=false;

			while(stop==false){

				double r1 = double(rand())/RAND_MAX; // random sample for p[k]
				double r2 = double(rand())/RAND_MAX; // random number for accepting
	
				double rmax = (c1[k]-1)/(c1[k]-1+c1[K-1]-1);
				double frac = pow(r1/rmax,c1[k]-1)*pow((1-r1)/(1-rmax),c1[K-1]-1);

				if(frac>r2){
					stop=true;
					r1=r1*(1.-pp+p[k]);
					p[K-1]+=p[k]-r1;
					pp-=p[k]-r1;  
					p[k]=r1;
					if(it%sampling_rate==0  &&  it>= Niter/10*9) sampling.push_back(r1);
				}
            } 
		}                                                                                                                                               		
		if(it%sampling_rate==0  &&  it>= Niter/10*9) sampling.push_back(1-pp);



              /* SAMPLING zeta ********************************************************************************************************/


                stop=false;

		double maxval=-100000;		
		double mx=0;

		int l=0;		

		while(l<10 && use_Potts==true && estimate_zeta==true) {
		
			double zeta1=0.5+0.05*l;

			double ZZ[K];
                        for(int k=0; k< K; k++) ZZ[k]=Zpart(N, NN[k], q, zeta1);

			double h=0;
                        for(int k=0; k< K; k++) h=h+NN[k]*log(ZZ[k]);

                        double val=(f1[0]-1)*log(zeta1)+(f1[1]-1)*log(1-zeta1)-h;

			if(val > maxval) { 
				maxval=val;
				mx=zeta1;
			}
			l++;
		}

	
                while(stop==false && use_Potts==true && estimate_zeta==true){


                        double r1 = double(rand())/RAND_MAX; // random sample for zeta
                        double r2 = double(rand())/RAND_MAX; // random number for accepting

                        double ZZ[K];
                        for(int k=0; k< K; k++) ZZ[k]=Zpart(N, NN[k], q, r1);
                        
			double h=0;
			for(int k=0; k< K; k++) h=h+NN[k]*log(ZZ[k]);

                        double val=(f1[0]-1)*log(r1)+(f1[1]-1)*log(1-r1)-h;


                        double frac = val - maxval;

 
			frac=exp(frac);

                        if(frac>r2){
                                stop=true;
                                if(it > 0) zeta = r1;
                                //if(it%sampling_rate==0) sampling.push_back(r1);
                        }

                }

		//if(use_Potts==true && it%sampling_rate==0  &&  it>= Niter/10*9) sampling.push_back(zeta);
 

		/****** SAMPLING Z *******************************************************************************************************/


		for(int i=0; i<N; i++){

			//if(i%1==0) cout << i << endl;
			if(fixed_Z==true) break;
       
			if(abs(zeta-1)<1E-5) {
				if(it%sampling_rate==0  &&  it>= Niter/10*9) sampling.push_back(Z[i]); 
				continue;
			}

			stop=false;

			double prob[K];
			double gg[K];
			double norm=0;
			double gmax=0;

			for(int k1=0; k1<K;k1++)	{	

				double g=0;

				if(use_Potts==true) {
			      		double n_in=0;
					for(int j=0; j<q;j++) {
						int index=Iin.at(q*i+j);
						if(Z[index]==k1) n_in=n_in+1.;
	      				}
					double m_in=0;
					for(int j=0; j<Iout_count[i];j++) {
						int index=Iout.at(Iout_track[i]+j); //Iout.at(N*i+j);
						if(index > -1 && Z[index]==k1) m_in=m_in+1.;
	      				}

					//cout << "nin " << n_in << " " << "min " << m_in << "nn " << Iout_count[i] << endl;

	                		//g=pow(zeta/(1-zeta),n_in+m_in)/Zpart(N,NN[k1],zeta,q);
					g=(n_in+m_in)*log(zeta/(1-zeta))-log(Zpart(N,NN[k1],zeta,q));					
					
					//g=g*pow(Zpart(N,NN[k1]-1,zeta,q)/Zpart(N,NN[k1],zeta,q),NN[k1]-1);
					g=g+log(Zpart(N,NN[k1]-1,zeta,q)/Zpart(N,NN[k1],zeta,q))*(NN[k1]-1);

				}
	
				if(g > gmax) gmax=g;

				gg[k1]=g;
				//cout << g << endl;

				//prob[k1]=p[k1]*d[k1]*pow(MU[i],-(d[k1]+1))*g;
				//prob[k1]=log(p[k1]*d[k1])-(d[k1]+1)*log(MU[i]);

				//norm+=prob[k1];

          		}

			for(int k1=0; k1<K;k1++)    gg[k1]=exp(gg[k1]-gmax); 

			for(int k1=0; k1<K;k1++)  {
				  prob[k1]=p[k1]*d[k1]*pow(MU[i],-(d[k1]+1))*gg[k1];
				  norm+=prob[k1];
			}

			for(int k1=0; k1<K;k1++)  prob[k1]=prob[k1]/norm;	

			while(stop==false){

				int r1 = rand()%K;
				double r2 = double(rand())/RAND_MAX;
				
				if(prob[r1] > r2) {
					stop=true;
					if(it%sampling_rate==0  &&  it>= Niter/10*9) sampling.push_back(r1);
					NN[Z[i]]-=1;
					a1[Z[i]]-=1;
					c1[Z[i]]-=1;
					V[Z[i]]-=log(MU[i]);
					b1[Z[i]]-=log(MU[i]);
					Z[i]=r1;
					NN[Z[i]]+=1;
					a1[Z[i]]+=1; 
					c1[Z[i]]+=1;
					V[Z[i]]+=log(MU[i]);
					b1[Z[i]]+=log(MU[i]);
					//cout << r1 << " ";
				}

			}

		}

		//cout << endl;


	/****** updating prior on zeta *******************************************************************************************************/

                N_in=0;
                for(int i=0; i< N; i++) {
                        int k=Z[i];
                        for(int j=0; j<q;j++) {
                                int index=Iin.at(q*i+j);
                                if(Z[index]==k) N_in=N_in+1.;
                    }
                }

		f1[0]=f[0]+N_in;
                f1[1]=f[1]+N*q-N_in;

/* likelihood ****************************************************************************************************/

		double lik0=0;
                for(int i=0; i< N; i++) lik0=lik0+log(p[Z[i]])+log(d[Z[i]])-(d[Z[i]]+1)*log(MU.at(i));
                
                double lik1=lik0+log(zeta/(1-zeta))*N_in;

                for(int k1=0; k1<K; k1++) lik1=lik1-NN[k1]*log(Zpart(N, NN[k1], zeta, q));
		



	/* save data **********************************************************************************************************/
	

		if(it%sampling_rate==0  &&  it>= Niter/10*9) sampling.push_back(lik0);
		if(it%sampling_rate==0  &&  it>= Niter/10*9) sampling.push_back(lik1);




	}


}
