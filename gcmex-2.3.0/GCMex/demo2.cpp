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
    mxArray *image_in_m, *mask_in_m, *c_out_m;
    const mwSize *dims;
    double *image, *mask, *c, sigma;
    int MXS, MYS, MZS, numdims;
    int x,y,z;

    /* Check for proper number of input and output arguments */
    if (nrhs != 3) {
        mexErrMsgIdAndTxt( "MATLAB:demo2:invalidNumInputs",
                "Wrong number of arguments\n\n\t Data usage : labelfield = demo2(image,mask,sigma).");
    }
    
    /* Check data type of input argument  */
    if (!(mxIsDouble(prhs[0])) || !(mxIsDouble(prhs[1])) || !(mxIsDouble(prhs[2]))){
        mexErrMsgIdAndTxt( "MATLAB:demo2:inputNotDouble",
                " Input arguments must be of type double.");
    }
    
//associate inputs
    image_in_m = mxDuplicateArray(prhs[0]);
    mask_in_m = mxDuplicateArray(prhs[1]);
    //d_in_m = mxDuplicateArray(prhs[2]); 
//figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
  MZS = (int)dims[0];  MYS = (int)dims[1];  MXS = (int)dims[2];
  // MZS = (int)dims[2];  MYS = (int)dims[0];  MXS = (int)dims[1];

//associate outputs
    c_out_m = plhs[0] = mxCreateDoubleMatrix(MYS*MXS* MZS,1,mxREAL);
    
//associate pointers
    image = mxGetPr(image_in_m);
    mask = mxGetPr(mask_in_m);
    sigma = mxGetScalar(prhs[2]);
    c = mxGetPr(c_out_m);

 int nbLabels = 3;   
///*do something*/
      GCoptimizationGeneralGraph gc(MXS*MYS*MZS,nbLabels);
  /* */
      
   for(z=0;z<MZS;z++)
    {  
   for(y=0;y<MYS;y++)
    {
          
         
       for(x=0;x<MXS;x++)
          {
              unsigned int pixelIndex =z*(MXS*MYS)+ y*MXS+x;
               //unsigned int pixelIndex = y*MXS+x;
              
              if ((int)(mask[pixelIndex])==0 ){
                
  
                  gc.setDataCost(pixelIndex,0,0.0);
                  gc.setDataCost(pixelIndex,1,10000.0);
                  gc.setDataCost(pixelIndex,2,10000.0);  
              }
              else if ((int)(mask[pixelIndex])== 59 ){
                  gc.setDataCost(pixelIndex,0,10000.0);
                  gc.setDataCost(pixelIndex,1,0.0);
                  gc.setDataCost(pixelIndex,2,10000.0);
              }else if ((int)(mask[pixelIndex])==195){
                  gc.setDataCost(pixelIndex,0,10000.0); 
                  gc.setDataCost(pixelIndex,1,10000.0);
                  gc.setDataCost(pixelIndex,2,0.0);  
              }
          }
    }
   }
     /* Connect horizontal neighboring pixels */
      /*
  for (int y = 0; y < MYS; y++ ){
    for (int  x = 1; x < MXS; x++ ){

        int pixelIndex = y*MXS+x;
        int pixelIndex2 = y*MXS+x-1;
        double gradient = (double)(fabs(image[pixelIndex]-image[pixelIndex2]));
        gc.setNeighbors(pixelIndex,pixelIndex2,(double)(1000*exp(-gradient/sigma))); 
    }
  }
    
  // Connect vertical neighboring pixels 
  for (int y = 1; y < MYS; y++ ){
    for (int  x = 0; x < MXS; x++){
        int pixelIndex = y*MXS+x;
        int pixelIndex2 = (y-1)*MXS+x;
        double gradient = (double)(fabs(image[pixelIndex]-image[pixelIndex2]));
        gc.setNeighbors(pixelIndex,pixelIndex2,(double)(1000*exp(-gradient/sigma))); 
        
    }
  }

*/
    
     /* Connect neighboring pixels */
    for(z=1;z<MZS;z++){
        for (int y = 1; y < MYS; y++ ){
            for (int  x = 1; x < MXS; x++ ){
               unsigned int pixelIndex1 = z*MXS*MYS+y*MXS+x;
                unsigned int pixelIndex2 = z*MXS*MYS+y*MXS+x-1;
                unsigned int pixelIndex3 = z*MXS*MYS+(y-1)*MXS+x;
                unsigned int pixelIndex4 = (z-1)*MXS*MYS+y*MXS+x;   
                //unsigned int pixelIndex1 = y*MXS+x;
               // unsigned int pixelIndex2 = y*MXS+x-1;
               // unsigned int pixelIndex3 = (y-1)*MXS+x;
             
                double gradient2 = (double)(fabs(image[pixelIndex1]-image[pixelIndex2]));
                double gradient3 = (double)(fabs(image[pixelIndex1]-image[pixelIndex3]));
                 double gradient4 = (double)(fabs(image[pixelIndex1]-image[pixelIndex4]));
                 
                gc.setNeighbors(pixelIndex1,pixelIndex2,(double)(1000*exp(-gradient2/sigma)));
                gc.setNeighbors(pixelIndex1,pixelIndex3,(double)(1000*exp(-gradient3/sigma)));
               gc.setNeighbors(pixelIndex1,pixelIndex4,(double)(1000*exp(-gradient4/sigma)));
                
            }
        }
    }
      
      
 /* go graph cut! */
 gc.expansion(9);

 /* save results */
  int i=0;
    for( int z=0;z<MZS;z++) {
        for( int y=0;y<MYS;y++){
            for( int x=0;x<MXS;x++,i++){
                c[i]=gc.whatLabel(i);
            }
        }
   }
}
