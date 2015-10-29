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
    mxArray *image_in_m, *mask_in_m, *c_out_m, *matrix_in_m;
    const mwSize *dims;
    const mwSize *dims_space, *dims_matrix;
    double *image, *mask, *c, sigma, *matrix;
    int MXS, MYS, MZS, numdims, Num_features, Num_points,Nb_points, nbLabels;
    int x,y,z;

    /* Check for proper number of input and output arguments */
    if (nrhs != 4) {
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
   matrix_in_m = mxDuplicateArray(prhs[2]);
   sigma = mxGetScalar(prhs[3]);
   
    //d_in_m = mxDuplicateArray(prhs[2]); 
//figure out dimensions
    dims = mxGetDimensions(prhs[1]);
    numdims = mxGetNumberOfDimensions(prhs[1]);
  MZS = (int)dims[2];  MYS = (int)dims[0];  MXS = (int)dims[1];
  // MZS = (int)dims[2];  MYS = (int)dims[0];  MXS = (int)dims[1];
    dims_matrix = mxGetDimensions(prhs[2]);
    Nb_points = (int)dims_matrix[1],   nbLabels = (int)dims_matrix[0];
  // assign number of features and number of points
  
  dims_space = mxGetDimensions(prhs[0]);
  Num_features = (int)dims_space[0];
  Num_points = (int)dims_space[1];
 // Num_features = 1;
  
//associate outputs
    c_out_m = plhs[0] = mxCreateDoubleMatrix(MYS*MXS* MZS,1,mxREAL);
    
//associate pointers
    image = mxGetPr(image_in_m);
    mask = mxGetPr(mask_in_m);
    matrix = mxGetPr(matrix_in_m);
   
    c = mxGetPr(c_out_m);

   
///*do something*/
      GCoptimizationGeneralGraph gc(MXS*MYS*MZS,nbLabels);
   
      
   for(z=0;z<MZS;z++)
    {  
  for(x=0;x<MXS;x++){
      for(y=0;y<MYS;y++)
    {
              unsigned int pixelIndex =z*(MXS*MYS)+ x*MYS+y;
               //unsigned int pixelIndex = y*MXS+x;
              
              if ((int)(mask[pixelIndex])==10 ){
                  gc.setDataCost(pixelIndex,0,0.0);
 
                  gc.setDataCost(pixelIndex,2,1000000.0);
 
              }else if ((int)(mask[pixelIndex])==2){
                  gc.setDataCost(pixelIndex,0,1000000.0); 

                  gc.setDataCost(pixelIndex,2,0.0);

              } 
              else
              {
                 for(int label=0;label<nbLabels;label++)
                  {
                        double datacost = matrix[pixelIndex * nbLabels + label];
                        gc.setDataCost(pixelIndex,label, datacost);
                  }

             }  
                
          }
    }
   }
   
    
    //Connect neighboring pixels 
for(z=1;z<MZS;z++){
   for(x=1;x<MXS;x++){
      for(y=1;y<MYS;y++){
            
               unsigned int pixelIndex1 = z*MXS*MYS+x*MYS+y;
                unsigned int pixelIndex2 = z*MXS*MYS+x*MYS+ y-1;
                unsigned int pixelIndex3 = z*MXS*MYS+ (x-1)*MYS+y;
                unsigned int pixelIndex4 = (z-1)*MXS*MYS+x*MYS+y;   
                //unsigned int pixelIndex1 = y*MXS+x;
               // unsigned int pixelIndex2 = y*MXS+x-1;
               // unsigned int pixelIndex3 = (y-1)*MXS+x;
             
                double gradient2 = 0;
                double gradient3 = 0;
                double gradient4 = 0;
                for (int f=0 ; f<Num_features;f++)
                {
                    unsigned int featureIndex1 = pixelIndex1 * Num_features + f;
                    unsigned int featureIndex2 = pixelIndex2 * Num_features + f;
                    unsigned int featureIndex3 = pixelIndex3 * Num_features + f;
                    unsigned int featureIndex4 = pixelIndex4 * Num_features + f;
                    gradient2 += (double)pow(image[featureIndex1] - image[featureIndex2],2);
                    gradient3 += (double)pow(image[featureIndex1] - image[featureIndex3],2);
                    gradient4 += (double)pow(image[featureIndex1] - image[featureIndex4],2);
                }
                    gradient2 = sqrt(gradient2);
                    gradient3 = sqrt(gradient3);   
                    gradient4 = sqrt(gradient4);     
                  //  char buffer [33];
                   // if (z==40 && y==40 && x==40)
                   // {
                   // itoa (gradient2,buffer,10);
                   // mexPrintf("\n");
                   // mexPrintf(buffer);
                   // itoa (gradient3,buffer,10);
                   // mexPrintf("\n");
                    //mexPrintf(buffer);
                   // itoa (gradient4,buffer,10);
                  //  mexPrintf("\n");
                  //  mexPrintf(buffer);
                  //  //}
                            
                 
                gc.setNeighbors(pixelIndex1,pixelIndex2,(double)(1000* exp(-gradient2/sigma)));
                gc.setNeighbors(pixelIndex1,pixelIndex3,(double)(1000* exp(-gradient3/sigma)));
                gc.setNeighbors(pixelIndex1,pixelIndex4,(double)(1000* exp(-gradient4/sigma)));
                
            }
        }
    }
      
      
// go graph cut!
 gc.expansion(9);

 // save results
  int i=0;
    for( int z=0;z<MZS;z++) {
       for(x=0;x<MXS;x++){
        for(y=0;y<MYS;y++, i++){
    
                c[i]=gc.whatLabel(i);
            }
        }
   }
  

}
