/* This is my HOG extractor writen based on the Matlab code. 
I do not perform any parameter checking, so do it yourself. 
The function assumes (check this or makre sure it is like this): 
	+ that the input image is of type double (no uint8)
	+ the blocks and cells are square
	+ the block and cell size is even
	+ there is no overlap between blocks
	*/

/*
Author: Vitomir Struc
Date: 31.8.2015 
*/

#include "mex.h" // tole ocitno more bit
#include <math.h>
#include <string.h>
#include <stdint.h>


/* globals */
int m = -1; // width of image - to je stevilo vrstic in je visina, spodaj je sirina
int n = -1; // height of image
int c = -1; // color information
uint8_t* Img = NULL;	// input image
double* pts = NULL;	// points of interest
double* features = NULL;	// output feature matrix
double* gaussian = NULL; // gausiian weights
double* x1y1 = NULL;
double* x1y2 = NULL;
double* x2y1 = NULL;
double* x2y2 = NULL;
double gMag = -1;
double gDir = -1;
double gDirBin = -1;
int gDirBinCell = -1;
int* cellX = NULL;
double* hist = NULL;
int pm = -1; // width of points
int pn = -1; // height of points

/* interpolation */
double x1y1z1 =  -1;
double x1y1z2 =  -1;
double x2y1z1 =  -1;
double x2y1z2 =  -1;
double x1y2z1 =  -1;
double x1y2z2 =  -1;
double x2y2z1 =  -1;
double x2y2z2 = -1;



/* some data i need */
int blocksizeInPixels = -1;
int featureSize = -1;
int halfSize = -1;
int numPoints = -1;
int validPointcount = -1;
int* validPointIdx = NULL;
double histRange = -1;

/* global hog settings */
int cell_size = -1;
int block_size = -1;
int overlap = -1;
int nbins = -1;
bool signedhist = false;

/* pixel access - I will only support grey-scale HOG, no color */
#define access_b(M,a,b) M[(a)+blocksizeInPixels*(b)]
#define PI 3.14159265
#define EPS 1.1921e-7

/* konstante */
double PI180 = 180/PI;
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


/* here we compute the gaussain weights */
void gaussianWeights()
{
	double sigma  = 0.5*blocksizeInPixels;	
	double siz = (blocksizeInPixels-1)/2.;
	double scale = 1./(2*sigma*sigma);
	double suma = 0; 
	double yc = 0;
	double xc = 0;
	double arg = 0;
	double maxv =  -100;
	int tmp =0;
	for(int i=0;i<blocksizeInPixels;++i)
	{
		yc = -siz+i;
		tmp = i*blocksizeInPixels;
		for(int j=0;j<blocksizeInPixels;++j)
		{
			xc = -siz+j;
			arg = -(xc*xc+yc*yc)*scale;
			gaussian[tmp+j] = exp(arg);
			suma += gaussian[tmp+j];
		}
	}

	if(suma !=0)
	{
		for(int i=0;i<blocksizeInPixels;++i)
		{
			tmp = i*blocksizeInPixels;
			for(int j=0;j<blocksizeInPixels;++j)
			{
				gaussian[tmp+j] /= suma;
			}
		}
	
	}

}


/* find bin center again */
void computeLowerHistBin1(double& x, double& binWidth, double& x1, double& X1)
{
	double invWidth = (double)(1./binWidth);
	int bin = (int)floor(x*invWidth-0.5);
	x1 = binWidth*(bin+0.5);
	X1 = (bin+2);
}


/* find bin center */
void computeLowerHistBin(double x, double& x1, int& X1)
{
	double invWidth = (double)(1./cell_size);
	int bin = (int)floor(x*invWidth-0.5);
	x1 = cell_size*(bin+0.5);
	X1 = (int)(bin+2);
}

/* precompute spatial histogram weights */
void spatialHistWeights()
{
	int width = cell_size*block_size;
	int height = width;
	double tmp = 0;
	int tmp1 = 0;
	double* xb = (double*) malloc( width*sizeof(double) );
	//memset(xb,0,width*sizeof(double));
	for(int i=0; i<width-0.5; ++i)// I will only have square cells and blocks - no need for cellY, yb
	{		
		computeLowerHistBin((double)(i+0.5), tmp, tmp1);
		cellX[i] = tmp1;
		xb[i] = 1-((double)(i+0.5)-tmp)/cell_size;	//to je zej wx1, wy1
	}

	for(int i=0; i<width-0.5; ++i)// I will only have square cells and blocks - no need for cellY, yb
	{	
		tmp1 = i*width;
		for(int j=0; j<height-0.5; ++j)
		{
			x1y1[tmp1+j] = xb[i]*xb[j];
			x1y2[tmp1+j] = xb[i]*(1-xb[j]);
			x2y1[tmp1+j] =(1- xb[i])*xb[j];
			x2y2[tmp1+j] = (1-xb[i])*(1-xb[j]);
		}

	}
	free(xb);
}

/* correct orientation scale */
inline void adjust_scale(double& com, double& hrange)
{
	if(com < 0)
		com += hrange;
}


/* compute gradients */
void hogGradient(int* roi)
{
	// for the gradients
	int gx = 0;
	int gy = 0;
	double pix1 = 0;
	double pix2 = 0;
	int roi_height = roi[1]+roi[3]; // roi is ok
	int roi_weight = roi[0]+roi[2];

	// for the range
	if(signedhist)
		histRange = 360;
	else
		histRange = 180;

	double binWidth = histRange/nbins; // tole je zdej ok
	int run_ind = 0;
	int cx = 0;
	int cy = 0;
	int z = 0;
	int kk =0;
	int kkm=0;
	int kkp=0;
	int hist_idx = 0;
	int tm=0;
	int tm1=0;
	double itm = 0;

	/* iterate through ROI */		
	for(int i=roi[0],ik=0; i<roi_weight; ++i,++ik)
	{
		cx = cellX[ik]-1;
		tm = (block_size+2)*(nbins+2)*cx;
		tm1 = (block_size+2)*(nbins+2)*(cx+1);
		kk = i*m;
		kkm = (i-1)*m;
		kkp = (i+1)*m;
		for(int j=roi[1],jk=0; j<roi_height; ++j,++jk)
		{
			cy = cellX[jk]-1;

			// x direction
			/*pix1 =  Img[(i-1)*m+j];
			pix2 =  Img[(i+1)*m+j];
			gx = (pix2-pix1);*/
			gx=(Img[kkp+j]-Img[kkm+j]);

			// y direction
			/*pix1 =  Img[(i)*m+j-1];
			pix2 =  Img[(i)*m+j+1];
			gy = -pix2+pix1;*/
			gy = -Img[kk+j+1]+Img[kk+j-1];

			pix1 = atan2((double)gy,(double)gx)*PI180;
			adjust_scale(pix1, histRange);
			run_ind = ik*blocksizeInPixels+(jk);
			gMag = sqrt((double)(gx*gx+gy*gy))*gaussian[run_ind];
			gDir = pix1;				
			computeLowerHistBin1(pix1, binWidth, pix2, itm);
			gDirBin = 1-(pix1-pix2)/binWidth;//to je wDir
			gDirBinCell = itm;

			//trilinear interpolation
			x1y1z1 = gDirBin*x1y1[run_ind]*gMag;
			x1y1z2 = (x1y1[run_ind]*gMag - x1y1z1);
			x2y1z1 = gDirBin*x2y1[run_ind]*gMag;
			x2y1z2 = (x2y1[run_ind]*gMag - x2y1z1);
			x1y2z1 = gDirBin*x1y2[run_ind]*gMag;
			x1y2z2 = (x1y2[run_ind]*gMag - x1y2z1);
			x2y2z1 = gDirBin*x2y2[run_ind]*gMag;
			x2y2z2 = (x2y2[run_ind]*gMag - x2y2z1);

			// do the binning
			z = gDirBinCell-1;
			hist_idx =  (tm+cy*(nbins+2)+z); // to je vse ok, tudi za 16x16 celice in 4x4 bloke
			hist[hist_idx] +=x1y1z1;

			hist_idx =  (tm+cy*(nbins+2)+z+1);
			hist[hist_idx]+=x1y1z2;

			hist_idx =  (tm+(cy+1)*(nbins+2)+z);
			hist[hist_idx] +=x1y2z1;

			hist_idx =  (tm+(cy+1)*(nbins+2)+z+1);
			hist[hist_idx] +=x1y2z2;

			hist_idx =  (tm1+(cy)*(nbins+2)+z);
			hist[hist_idx]+=x2y1z1;

			hist_idx =  (tm1+(cy)*(nbins+2)+z+1);
			hist[hist_idx]+=x2y1z2;
				
			hist_idx =  (tm1+(cy+1)*(nbins+2)+z);
			hist[hist_idx] +=x2y2z1;
				
			hist_idx =  (tm1+(cy+1)*(nbins+2)+z+1);
			hist[hist_idx] +=x2y2z2;

		}
	}

	// wrap it
	int hist_idx1 = 0;
	for(int i=0; i<block_size+2; ++i)
	{
			cx = i;
			tm = (block_size+2)*(cx)*(nbins+2);
			for(int j=0; j<block_size+2; ++j)
			{
				cy = j;
				hist_idx =  (tm+(cy)*(nbins+2)+1);
				hist_idx1 =  (tm+(cy)*(nbins+2)+nbins+1);
				hist[hist_idx] += hist[hist_idx1];

				hist_idx =  (tm+(cy)*(nbins+2)+nbins);
				hist_idx1 =  (tm+(cy)*(nbins+2)+0);
				hist[hist_idx] += hist[hist_idx1];				
			}
	}
}

/* Compute final HOG */
void extractHog(int* roi, int pointNo)
{
	int r = (int)blocksizeInPixels;
	int c = (int)blocksizeInPixels;
	int nCells = (int)block_size;
	int numCellsPerBlock = nCells*nCells;
	int roi_height = roi[1]+roi[3]; // roi is ok
	int roi_weight = roi[0]+roi[2];
	
	//we have the histogram - wrap it
	int z2 = 1;
	int run_ind = (roi_weight-1-roi[0])*blocksizeInPixels+(roi_height-1-roi[1]);
	int hist_idx=0;
	int hist_idx1=0;
	int cx=0;
	int cy =0;
	int tm = 0;
	int tm1 = 0;

	////normalize it
	double sum =0;
	double sum1=0;
	for(int i=1; i<block_size+1; ++i)
	{
			cx = i;
			tm=(block_size+2)*(cx)*(nbins+2);
			for(int j=1; j<block_size+1; ++j)
			{
				cy = j;
				tm1 = tm+(cy)*(nbins+2);
				for(int k=1;k<nbins+1;++k)
				{
					hist_idx =  (tm1+k);
					sum += (hist[hist_idx]*hist[hist_idx]);
				}
			}
	}
	for(int i=1; i<block_size+1; ++i)
	{
			cx = i;
			tm=(block_size+2)*(cx)*(nbins+2);
			for(int j=1; j<block_size+1; ++j)
			{
				cy = j;
				tm1 = tm+(cy)*(nbins+2);
				for(int k=1;k<nbins+1;++k)
				{
					hist_idx =  (tm1+k);
					hist[hist_idx]/=(sqrt(sum)+EPS);
					if(hist[hist_idx]>0.2)
						hist[hist_idx] = 0.2;
					sum1+=(hist[hist_idx]*hist[hist_idx]);
				}
			}
	}
	int tm2=0;
	int tm3=0;
	for(int i=1; i<block_size+1; ++i)
	{
			cx = i;
			tm=(block_size+2)*(cx)*(nbins+2);
			tm1=(block_size)*(cx-1)*(nbins);
			for(int j=1; j<block_size+1; ++j)
			{
				cy = j;
				tm2=tm+(cy)*(nbins+2);
				tm3=tm1+(cy-1)*(nbins);
				for(int k=1;k<nbins+1;++k)
				{
					hist_idx =  (tm2+k);
					hist[hist_idx]/=(sqrt(sum1)+EPS);
					z2 =  pointNo+(numPoints)*(tm3+k-1);
					features[z2]=hist[hist_idx];
				}
			}
	}



}

/* main HOG function */
void  extractHOGFromPoints()
{
	/* lets precompute some stuff */
	blocksizeInPixels = cell_size*block_size;
	featureSize = block_size*block_size*nbins;
	halfSize = blocksizeInPixels/2; // check that blocksizeInPixels is even
    validPointcount = 0;
	int roi[4] = {1,1,blocksizeInPixels,blocksizeInPixels};

	for(int i=0; i<(int)numPoints; ++i)
	{
		roi[0] = pts[0*(int)numPoints+i]-1-halfSize;
		roi[1] = pts[1*(int)numPoints+i]-1-halfSize;

		//process if ROI is fully in image
		if(roi[0]>=1 && roi[1]>=1 && roi[0]+roi[2]<n && roi[1]+roi[3]<m)//we only execute the function when all points are in the image - so I can quickly pad it otherwise
		{
			validPointcount++;
		}
	}

	if(validPointcount == numPoints)
	{
		/* Create gaussian weights */
		gaussian = (double*) malloc( blocksizeInPixels*blocksizeInPixels*sizeof(double) );
		//memset(gaussian,0,blocksizeInPixels*blocksizeInPixels*sizeof(double));
		gaussianWeights(); // here we set "gaussian" 

		/* Precompute spatial histogram weights */
		x1y1 = (double*) malloc( blocksizeInPixels*blocksizeInPixels*sizeof(double) );
		//memset(x1y1,0,blocksizeInPixels*blocksizeInPixels*sizeof(double));
		x2y1 = (double*) malloc( blocksizeInPixels*blocksizeInPixels*sizeof(double) );
		//memset(x2y1,0,blocksizeInPixels*blocksizeInPixels*sizeof(double));
		x1y2 = (double*) malloc( blocksizeInPixels*blocksizeInPixels*sizeof(double) );
		//memset(x1y2,0,blocksizeInPixels*blocksizeInPixels*sizeof(double));
		x2y2 = (double*) malloc( blocksizeInPixels*blocksizeInPixels*sizeof(double) );
		//memset(x2y2,0,blocksizeInPixels*blocksizeInPixels*sizeof(double));
		cellX = (int*) malloc( blocksizeInPixels*sizeof(int) );
		//memset(cellX,0,blocksizeInPixels*sizeof(int));
		spatialHistWeights(); // here we precompute the weights for the trilinear interpolation
	
		hist = (double*) malloc( (nbins+2)*(block_size+2)*(block_size+2)*sizeof(double) );
	

		for(int i=0; i<(int)numPoints; ++i)
		{
			roi[0] = pts[0*(int)numPoints+i]-1-halfSize;
			roi[1] = pts[1*(int)numPoints+i]-1-halfSize;
			memset(hist,0,(nbins+2)*(block_size+2)*(block_size+2)*sizeof(double));		
			hogGradient(roi);
			extractHog(roi,i);		
		}
		free(gaussian);
		free(x1y1);
		free(x2y1);
		free(x1y2);
		free(x2y2);
		free(cellX);
		free(hist);
	}
}


//to je zej vstopna tocka 
void mexFunction(int nlhs, mxArray *plhs[],  /* Output variables */
			int nrhs, const mxArray *prhs[]) /* Input variables */
{
	// preverimo samo, da imam dovolj parametrov
	if( nrhs!=7 ) 
		mexErrMsgTxt("Always 7 input arguments required."); 

	/* Input 1 - the image */
	get_dimensions( prhs[0], m,n,c);
	Img = (unsigned char *)mxGetPr(prhs[0]);

	/* Input 2 - interest points */
	int s = 0;
	get_dimensions( prhs[1], numPoints,pn,s);//ce damo not 68x2 tocke, dobimo pm=68, pn=2
	pts = mxGetPr(prhs[1]);

	/* Input 3 - cell size */
	cell_size = (int) *mxGetPr(prhs[2]);

	/* Input 4 - block size */
	block_size = (int) *mxGetPr(prhs[3]);

	/* Input 5 - overlap */
	overlap = (int) *mxGetPr(prhs[4]);

	/* Input 6 - overlap */
	nbins = (int) *mxGetPr(prhs[5]);

	/* Input 7 - signed hist */
	signedhist = *mxGetPr(prhs[6])>0.5;

	/* Output 1 - HOG s */
	if(nlhs==0 || nlhs==1 || nlhs ==2)
	{
		plhs[0] =  mxCreateDoubleMatrix(numPoints,block_size*block_size*nbins, mxREAL); 
		features = mxGetPr(plhs[0]); // I use this in all functions	

		/* compute descriptor in subroutine*/
		extractHOGFromPoints();

		if(nlhs==2)
		{
			plhs[1] =  mxCreateDoubleScalar(validPointcount);
		}
	}
	else 
		mexErrMsgTxt("Only two output arguments supported.");	
	return;
}