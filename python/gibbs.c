#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <sstream>
//#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
//#include <algorithm>
//#include <numeric>
//#include <assert.h>
#include <Python.h>
#include "numpy/arrayobject.h"
#include "gibbs.h"



/* #### Globals #################################### */

/* ==== Set up the methods table ====================== */
static PyMethodDef _gibbsMethods[] = {
	{"GibbsSampling", GibbsSampling, METH_VARARGS},
	{NULL, NULL}     /* Sentinel - marks the end of this structure */
};

 static struct PyModuleDef _gibbs =
{

   PyModuleDef_HEAD_INIT,
   "cModPyDem", /* name of module */
    "",          /* module documentation, may be NULL */
    -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
   _gibbsMethods
};


PyMODINIT_FUNC PyInit__gibbs(void)
{
    import_array();  // Must be present for NumPy.  Called first after above line.
    return PyModule_Create(&_gibbs);
}


/* ==== Create 1D Carray from PyArray ======================
    Assumes PyArray is contiguous in memory.             */
double *pyvector_to_Carrayptrs(PyArrayObject *arrayin)  {
	int i,n;
	
	n=arrayin->dimensions[0];
	return (double *) arrayin->data;  /* pointer to arrayin data as double */
}

/*
long int *pyvector_to_Carrayiptrs(PyArrayObject *arrayin)  {
	int i,n;
	
	n=arrayin->dimensions[0];
	return (long int *) arrayin->data;  // pointer to arrayin data as double 
}
*/

// binomial coefficient
double binom(int N, int q) {

 double s=1;
 int q1; 

 if(q==0) return s;
 
 else{ 
 	for(q1=0; q1<q; q1++) {
              s=s*(N-q1)/(q1+1);
   	}
	return s; 
 }

}


// partition function Z
double Zpart(int N, int N1, double zeta, int q){

        double s=0;
	int q1;

        for(q1=0; q1<=q; q1++) {
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





static PyObject *GibbsSampling(PyObject *self, PyObject *args) {

//void Gibbs_Potts(int Niter, int K, bool fixed_Z, bool use_Potts, bool estimate_zeta, int q, vector <double> MU, vector<int> Iin, vector<int> Iout, int Iout_count[], int Iout_track[],  double a[], double b[], double c[], double f[], double zeta, vector <double>& sampling, int sampling_rate, int replica){

  
        int Niter, K, fixed_Z, use_Potts, estimate_zeta, q, sampling_rate, replica; 
	double zeta, burn_in;

        PyArrayObject *vecMU, *vecIin, *vecIout, *vecIout_count, *vecIout_track, *veca, *vecb, *vecc, *vecf, *vecsampling; 

	double *MU;
	double *Iin;
	double *Iout;
	double *Iout_count;
	double *Iout_track;
	double *sampling; 		
	double *a, *b, *c,*f;


        /* Parse tuples separately since args will differ between C fcns */
        if (!PyArg_ParseTuple(args, "iiiiiididiO!O!O!O!O!O!O!O!O!O!", &Niter, &K, &fixed_Z, &use_Potts, &estimate_zeta, &q, &zeta, &sampling_rate, &burn_in, &replica, &PyArray_Type, &vecMU, &PyArray_Type, &vecIin,&PyArray_Type, &vecIout,&PyArray_Type, &vecIout_count,&PyArray_Type, &vecIout_track,&PyArray_Type, &veca, &PyArray_Type, &vecb, &PyArray_Type, &vecc, &PyArray_Type, &vecf, &PyArray_Type, &vecsampling))  return NULL;
	if (NULL == vecMU)  return NULL;	
	if (NULL == vecIin)  return NULL;
	if (NULL == vecIout)  return NULL;
	if (NULL == vecIout_count)  return NULL;
	if (NULL == vecIout_track)  return NULL;
	if (NULL == veca)  return NULL;
	if (NULL == vecb)  return NULL;
	if (NULL == vecc)  return NULL;
	if (NULL == vecf)  return NULL;
	if (NULL == vecsampling)  return NULL;


	MU=pyvector_to_Carrayptrs(vecMU);
	Iin=pyvector_to_Carrayptrs(vecIin);
	Iout=pyvector_to_Carrayptrs(vecIout);
	Iout_count=pyvector_to_Carrayptrs(vecIout_count);
	Iout_track=pyvector_to_Carrayptrs(vecIout_track);
	a=pyvector_to_Carrayptrs(veca);
	b=pyvector_to_Carrayptrs(vecb);
	c=pyvector_to_Carrayptrs(vecc);
	f=pyvector_to_Carrayptrs(vecf);
	sampling=pyvector_to_Carrayptrs(vecsampling);

	int N=vecMU->dimensions[0];

	srand(0);

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

	int i,j,k;

	for(k=0; k<K; k++) { 
		V[k]=0;  
		NN[k]=0;
		d[k]=1;
		//sampling[k]=d[k];
	}

	for(k=0; k<K; k++) {
		p[k]=1./K;
		//sampling[K+k]=p[k];
	}

	pp=(K-1.0)/(double)K;

 	
	for(i=0; i<N; i++){
		if(fixed_Z==0) { 
			int z=rand()%K;
			Z[i]=z;
		}
		else Z[i]=0;

		V[Z[i]]=V[Z[i]]+log(MU[i]);
		NN[Z[i]]+=1;
		//sampling[2*K+i]=Z[i];
	}

	for(k=0; k<K; k++) {
		a1[k]=a[k]+NN[k];
		b1[k]=b[k]+V[k];
		c1[k]=c[k]+NN[k];
	}

	double N_in=0;
	for(i=0; i< N; i++) {
		int k=Z[i];
		for(j=0; j<q;j++) {
			int index=(int)Iin[q*i+j];
			if(Z[index]==k) N_in=N_in++;
	    }			
	}	


	f1[0]=f[0]+N_in; 
	f1[1]=f[1]+N*q-N_in;

        //sampling[2*K+N]=0;
        //sampling[2*K+N+1]=0;

/******************* ITERATIONS **********************************************************************************************************/
	
	int it;
	int prev=0;
	int prev1=0;

	for(it=0; it<Niter; it++){

		int par=0;

		if(it%100==0) printf("it= %d \n",it);

		int stop;


		/* SAMPLING d***********************************************************************************************************/

		for(k=0; k<K; k++){
	        
			stop=0;


			while(stop==0){
	   
				double r1 = (double)rand()/RAND_MAX*200;
				double r2 = (double)rand()/RAND_MAX;
	      
				double rmax = (a1[k]-1)/b1[k];
				double frac;

				if(a1[k]-1>0) frac = exp(-b1[k]*(r1-rmax)-(a1[k]-1)*(log(rmax)-log(r1))); 
				else frac = exp(-b1[k]*r1);


				if(frac>r2){  
					stop=1;
					if(it%sampling_rate==0 &&  it>= Niter*burn_in) sampling[prev1+k]=r1;
					d[k]=r1; 
				}
			}
		}

		par+=K;

		/* SAMPLING p **********************************************************************************************************/


		for(k=0; k<K-1; k++) {

			stop=0;

			while(stop==0){

				double r1 = (double)rand()/RAND_MAX; // random sample for p[k]
				double r2 = (double)rand()/RAND_MAX; // random number for accepting
	
				double rmax = (c1[k]-1)/(c1[k]-1+c1[K-1]-1);
				double frac = pow(r1/rmax,c1[k]-1)*pow((1-r1)/(1-rmax),c1[K-1]-1);

				if(frac>r2){
					stop=1;
					r1=r1*(1.-pp+p[k]);
					p[K-1]+=p[k]-r1;
					pp-=p[k]-r1;  
					p[k]=r1;
					if(it%sampling_rate==0  &&  it>= Niter*burn_in) sampling[prev1+par+k]=r1;
				}
                        } 
		}                                                                                                                                               		
		if(it%sampling_rate==0  &&  it>= Niter*burn_in)  sampling[prev1+par+K-1]=1-pp;

		par+=K;


              /* SAMPLING zeta ********************************************************************************************************/


                stop=0;

		double maxval=-100000;		
		double mx=0;

		int l=0;		

		while(l<10 && use_Potts==1 && estimate_zeta==1) {

		
			double zeta1=0.5+0.05*l;

			double ZZ[K];
                        for(k=0; k< K; k++) ZZ[k]=Zpart(N, NN[k], q, zeta1);

			double h=0;
                        for(k=0; k< K; k++) h=h+NN[k]*log(ZZ[k]);

                        double val=(f1[0]-1)*log(zeta1)+(f1[1]-1)*log(1-zeta1)-h;

			if(val > maxval) { 
				maxval=val;
				mx=zeta1;
			}
			l++;
		}

	
                while(stop==0 && use_Potts==1 && estimate_zeta==1){



                        double r1 = (double)rand()/RAND_MAX; // random sample for zeta
                        double r2 = (double)rand()/RAND_MAX; // random number for accepting

			double ZZ[K];
                        for(k=0; k< K; k++) ZZ[k]=Zpart(N, NN[k], q, r1);
                        
			double h=0;
			for(k=0; k< K; k++) h=h+NN[k]*log(ZZ[k]);

                        double val=(f1[0]-1)*log(r1)+(f1[1]-1)*log(1-r1)-h;


                        double frac = val - maxval;

 
			frac=exp(frac);

                        if(frac>r2){
                                stop=1;
                                if(it > 0) zeta = r1;
				if(it%sampling_rate==0  &&  it>= Niter*burn_in) sampling[prev1+par]=r1;
                        }
			par+=1;
                }

		//if(use_Potts==true && it%sampling_rate==0  &&  it>= Niter/10*9) sampling.push_back(zeta);



/****** SAMPLING Z *******************************************************************************************************/


		for(i=0; i<N; i++){

			if(fixed_Z==1) break;

			if(sqrt((zeta-1)*(zeta-1))<1E-5) { 
				if(it%sampling_rate==0  &&  it>= Niter*burn_in) sampling[prev1+par+i]=(double)Z[i]; 
				continue;
			}

			stop=0;

			double prob[K];
			double gg[K];
			double norm=0;
			double gmax=0;

			int k1; 

			for(k1=0; k1<K;k1++)	{	

				double g=0;

				if(use_Potts==1) {

			      		double n_in=0;
					for(j=0; j<q;j++) {
						int index=(int)Iin[q*i+j];
						if(Z[index]==k1) n_in=n_in+1.;
	      				}
					double m_in=0;
					for(j=0; j<(int)Iout_count[i];j++) {
						int index=(int)Iout[(int)Iout_track[i]+j]; //Iout.at(N*i+j);
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

			

			for(k1=0; k1<K;k1++)    gg[k1]=exp(gg[k1]-gmax); 

			for(k1=0; k1<K;k1++)  {
				  prob[k1]=p[k1]*d[k1]*pow(MU[i],-(d[k1]+1))*gg[k1];
				  norm+=prob[k1];
			}

			for(k1=0; k1<K;k1++)  prob[k1]=prob[k1]/norm;	

			while(stop==0){

				int r1 = rand()%K;
				double r2 = (double)rand()/RAND_MAX;
				
				if(prob[r1] > r2) {
					stop=1;
					if(it%sampling_rate==0  &&  it>= Niter*burn_in) sampling[prev1+par+i]=(double)r1;
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
				}

			}

		}

		par+=N;
		
        /****** updating prior on zeta *******************************************************************************************************/

                N_in=0;
                for(i=0; i< N; i++) {
                        int k=Z[i];
                        for(j=0; j<q;j++) {
                                int index=(int)Iin[q*i+j];
                                if(Z[index]==k) N_in=N_in+1.;
                    }
                }

		f1[0]=f[0]+N_in;
                f1[1]=f[1]+N*q-N_in;

/* likelihood ****************************************************************************************************/

		double lik0=0;
                for(i=0; i< N; i++) lik0=lik0+log(p[Z[i]])+log(d[Z[i]])-(d[Z[i]]+1)*log(MU[i]);
                
                double lik1=lik0+log(zeta/(1-zeta))*N_in;

                for(k=0; k<K; k++) lik1=lik1-NN[k]*log(Zpart(N, NN[k], zeta, q));
		



	/* save data **********************************************************************************************************/
	

		if(it%sampling_rate==0  &&  it>= Niter*burn_in) sampling[prev1+par]=lik0;
		if(it%sampling_rate==0  &&  it>= Niter*burn_in) sampling[prev1+par+1]=lik1;

		if(it%sampling_rate==0 && it>= Niter*burn_in) prev++;

		prev1=prev*(par+2);


	}


	return Py_BuildValue("i", 1);


}
