#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
     
    int numPoints, numConfigs;
    double a11, a12, a13, a21, a22, a23, a31, a32, a33;
    
    // input variables
    double *configs;
    double *inputXs_src, *inputYs_src;
    double *inputXs_tgt, *inputYs_tgt;
    
      
    // output variables
    double *distances;
    
    
    //// checking inputs
    //if (nrhs != 5)
    //	mexErrMsgTxt("The number of input arguments must be 5.");
    //if (nlhs != 1)
    //	mexErrMsgTxt("The number of output arguments must be 5.");
    
    
    // inputs:     dimensions
    // ------	    -----------
    // configs      numConfigs x 8
    // inputXs_src           1 x numPoints
    // ys1			1 x numPoints
    
    // outputs:     dimensions
    // -------	    -----------
    // distances    numPoints x numConfigs
    
    
    /* Find the dimensions of the data */
    numConfigs = mxGetN(prhs[0]);
    numPoints = mxGetN(prhs[1]);       
    
    /* Create an mxArray for the output data */
    plhs[0] = mxCreateDoubleMatrix(numPoints, numConfigs, mxREAL );
    
    
    /* Retrieve the input data */
    configs = mxGetPr(prhs[0]);
    inputXs_src = (double*)mxGetPr(prhs[1]);
    inputYs_src = (double*)mxGetPr(prhs[2]);
    inputXs_tgt = (double*)mxGetPr(prhs[3]);
    inputYs_tgt = (double*)mxGetPr(prhs[4]);
   
        
    /* Create a pointer to the output data */
    distances = mxGetPr(plhs[0]);

        
    // MAIN LOOP

	  // helper variables
    double 	translatedPoint_x, translatedPoint_y;
	double 	tagetPoint_norm;
	double 	*ptrXsrc, *ptrYsrc;
	double 	*ptrXtgt, *ptrYtgt;
	double  *ptrDs;// = distances;
    int		i,j;
#pragma omp parallel for if (numConfigs > 1000) num_threads(omp_get_num_procs()) default(shared) private(j,i,ptrDs, translatedPoint_x, translatedPoint_y, ptrXsrc, ptrYsrc, ptrXtgt, ptrYtgt, tagetPoint_norm, a11,a12,a13,a21,a22,a23,a31,a32,a33)
    for (i = 0 ; i < numConfigs ; i++) {
        
//        if (i%100000==0)
//            mexPrintf("MAIN LOOP: config %d out of %d\n", i+1, numConfigs);
//		if (i==709416)
//			mexPrintf("MAIN LOOP: config %d out of %d\n",i+1,numConfigs);
        
        a11 = configs[9*i];
        a12 = configs[9*i+1];
        a13 = configs[9*i+2];
        a21 = configs[9*i+3];
        a22 = configs[9*i+4];
        a23 = configs[9*i+5];
        a31 = configs[9*i+6];
        a32 = configs[9*i+7];
        a33 = configs[9*i+8];

		ptrXsrc = inputXs_src;
        ptrYsrc = inputYs_src;
		ptrXtgt = inputXs_tgt;
        ptrYtgt = inputYs_tgt;
		ptrDs = distances + numPoints * i;
        for (j = 0; j < numPoints ; j++) 
		{
			tagetPoint_norm =  a31*(*ptrXsrc) + a32*(*ptrYsrc) + a33; // includes rounding
            translatedPoint_x = (a11*(*ptrXsrc) + a12*(*ptrYsrc) + a13) / tagetPoint_norm; 
            translatedPoint_y = (a21*(*ptrXsrc) + a22*(*ptrYsrc) + a23) / tagetPoint_norm; 

			*ptrDs = sqrt((*ptrXtgt-translatedPoint_x)*(*ptrXtgt-translatedPoint_x) + (*ptrYtgt-translatedPoint_y)*(*ptrYtgt-translatedPoint_y));

            //targetInd = max(0,(translatedPoint_y - 1)*w2 + translatedPoint_x - 1); // -1 is for c
            //*ptrDs = img2[targetInd];
                
            ptrXsrc++;
            ptrYsrc++;
            ptrXtgt++;
            ptrYtgt++;
			ptrDs++;
        }        

    }
    
    
}