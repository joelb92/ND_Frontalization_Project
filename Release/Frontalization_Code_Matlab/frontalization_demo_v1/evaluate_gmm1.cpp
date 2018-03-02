/*
Author: Vitomir Struc
Date: 31.8.2015 
*/

#include "mex.h" // tole ocitno more bit
#include <math.h>
#include <string.h>

int dim=-1;
int num_data=-1;
int c=-1;
int d=-1;
int e=-1;
int ncomp=-1;
double* X = NULL;
double* M = NULL;
double* C = NULL;
double* y = NULL;

/* macros */
#define PI 3.14159265
#define EPS 1.1921e-7

/* Read dimensions of array */
inline void get_dimensions( const mxArray* arr, int& a, int& b, int& c )
{
	int nd = mxGetNumberOfDimensions(arr);
	if( nd==3 )
	{
		a = mxGetDimensions(arr)[0];
		b = mxGetDimensions(arr)[1];
		c = mxGetDimensions(arr)[2];
	}
	else if( nd==2 )
	{
		a = mxGetM(arr);
		b = mxGetN(arr);
		c = 1;
	}
	else
		mexErrMsgTxt("Only 2D and 3D arrays supported."); 
}

/* helpers */
inline void display_message(const char* mess, int v)
{
	char str[128];
	sprintf(str, mess, v);
	mexWarnMsgTxt(str);
}

/* helpers */
inline void display_message_f(const char* mess, double v)
{
	char str[128];
	sprintf(str, mess, v);
	mexWarnMsgTxt(str);
}

void evaluate()
{
	double  con = con = dim*log(2*PI);
	double sum=0;
	double clogs = 0;
	double ysum = 0;
	int idx = 0;
	double maxy = -100000;
	//display_message("dim=%d", dim);
	for(int i=0;i<ncomp;++i)
	{
		sum=0;
		clogs = 0;
		idx = i*dim;
		for(int j=0;j<dim;++j)
		{
			sum   += (X[j]-M[idx+j])*(X[j]-M[idx+j])/(C[idx+j]+EPS);
			clogs += log(C[idx+j]+EPS);
		}
		y[i] = (-0.5*sum)-0.5*(con+clogs); 
		//display_message_f("y[j]=%.8f", y[i]);
		//display_message_f("clogs=%.8f", clogs);
		//display_message_f("sum=%.8f", sum);
		if(maxy<y[i])
			maxy = y[i]; 
	}
	
	for(int i=0;i<ncomp;++i)
	{
		clogs = y[i] - maxy;
		ysum += exp(clogs);		
	}
	ysum = maxy+log(ysum);
	
	for(int i=0;i<ncomp;++i)
	{
		y[i] = exp(y[i]-ysum);//+EPS;
		if(y[i]!=y[i]){
			y[i]=0;}
	}
}


//to je zej vstopna tocka 
void mexFunction(int nlhs, mxArray *plhs[],  /* Output variables */
			int nrhs, const mxArray *prhs[]) /* Input variables */
{
	

	/* Input 1 - the data */
	get_dimensions( prhs[0], dim,num_data,c);
	X = mxGetPr(prhs[0]);

	
	/* Input 2 - the mean */
	get_dimensions( prhs[1], d,ncomp,c);
	M = mxGetPr(prhs[1]);

	/* Input 3 - covariance */
	get_dimensions( prhs[2], d,e,c);
	C = mxGetPr(prhs[2]);

	/* Output 1 - y */
	plhs[0] =  mxCreateDoubleMatrix(ncomp,num_data, mxREAL); 
	y = mxGetPr(plhs[0]);

	/* compute distances in subroutine*/
	evaluate();

}