/*********************************************************************
 * Demo.cpp
 *
 * This file shows the basics of setting up a mex file to work with
 * Matlab.  This example shows how to use 2D matricies.  This may
 * 
 * Keep in mind:
 * <> Use 0-based indexing as always in C or C++
 * <> Indexing is column-based as in Matlab (not row-based as in C)
 * <> Use linear indexing.  [x*dimy+y] instead of [x][y]
 *
 * For more information, see my site: www.shawnlankton.com
 * by: Shawn Lankton
 *
 ********************************************************************/
#include <matrix.h>
#include "GCoptimization.h"
#include <mex.h>   
#include <stdlib.h>
#include <cmath>
#include <math.h>
/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//declare variables
    mxArray *image_in_m,  *c_out_m , *matrix_in_m;
    const mwSize *dims , *size;
    double *image, *c, *matrix, beta;
    int MXS, MYS, numdims;
    int x,y;
      int Nb_points, nbLabels ;

    /* Check for proper number of input and output arguments */
    if (nrhs != 3) {
        mexErrMsgIdAndTxt( "MATLAB: markov_network_multilable:invalidNumInputs",
               "Wrong number of arguments\n\n\t Data usage : labelfield = demo2(image,mask,sigma).");
    }
  
    /* Check data type of input argument  */
    if (!(mxIsDouble(prhs[0])) || !(mxIsDouble(prhs[1]))  || !(mxIsDouble(prhs[2]))){
       mexErrMsgIdAndTxt( "MATLAB: markov_network_multilable:inputNotDouble",
               " Input arguments must be of type double.");
  }
    

//associate inputs
 image_in_m = mxDuplicateArray(prhs[0]);
 
    //d_in_m = mxDuplicateArray(prhs[2]); 
//figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    MYS = (int)dims[0]; MXS = (int)dims[1];

//associate outputs
    c_out_m = plhs[0] = mxCreateDoubleMatrix(MYS,MXS,mxREAL);
    
    beta = mxGetScalar(prhs[2]);
//associate pointers
    image = mxGetPr(image_in_m);
  matrix_in_m = mxDuplicateArray(prhs[1]);
    matrix = mxGetPr(matrix_in_m);
   
    c = mxGetPr(c_out_m);

     size = mxGetDimensions(prhs[1]);
     Nb_points = (int)size[1],   nbLabels = (int)size[0];
    char buffer [33];
       
///*do something*/
    //int nbLabels = 256;
       GCoptimizationGeneralGraph gc(MXS*MYS,nbLabels);
   
        int ii;
    for(y=0;y<MYS;y++)
    {
          for(x=0;x<MXS;x++)
          {
            int pixelIndex = y*MXS+x;
            for(int label=0;label<nbLabels;label++)
            {
                int datacost = matrix[pixelIndex * nbLabels + label];
                gc.setDataCost(pixelIndex,label,datacost);
                /*if (pixelIndex == 0)
                    {
                    itoa (datacost,buffer,10);
                    mexPrintf("\n");
                    mexPrintf(buffer);
                    }*/
            }
               if (pixelIndex < 1000){
               itoa (image[pixelIndex],buffer,10);
                mexPrintf("  ");
                  mexPrintf(buffer); } 
          }
    }
        
        for ( int l1 = 0; l1 < nbLabels; l1++ )
            for (int l2 = 0; l2 < nbLabels; l2++ ){
                    float cost =  beta*abs(l1-l2);
                    gc.setSmoothCost(l1,l2,(int)cost);
            }     
        
         /* Connect horizontal neighboring pixels */
        for (int y = 0; y < MYS; y++ ){
                for (int  x = 1; x < MXS; x++ ){
                        int pixelIndex1 = y*MXS+x;
                        int pixelIndex2 = y*MXS+x-1;
                         gc.setNeighbors(x+y*MXS,x-1+y*MXS); 
                         
                }
        }

        /* Connect vertical neighboring pixels */
        for (int y = 1; y < MYS; y++ ){
                for (int  x = 0; x < MXS; x++){

                        int pixelIndex1 = y*MXS+x;
                        int pixelIndex2 = (y-1)*MXS+x; 
                        gc.setNeighbors(x+y*MXS,x+(y-1)*MXS);
                }
        }

       

  /* On lance l'algorithme de graph cut */
        gc.expansion(2);


 /* save results */
    int i=0;
    for( int y=0;y<MYS;y++){
        for( int x=0;x<MXS;x++,i++){
            c[i]=gc.whatLabel(i);
        }
    }
      
} 
