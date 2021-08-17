/* Change this to reflect the apropriate header file */
#include <math.h>
#include <limits.h>
#include <time.h>
#include "mex.h"
#include "matrix.h"

/*
 * composite_likelihood.c -
 * This is a helper function and is part of the MFE toolbox
 * You should be able to compile it on any platform.
 *
 * Author: Kevin Sheppard
 * kevin.sheppard@economics.ox.ac.uk
 * Revision: 1    Date: 4/17/2012
 */

void composite_likelihood_core(double *S, double *X, double *indices, size_t q, size_t m, size_t n, double *ll) {
   /* q = size(indices,1);
    * // [m,n] = size(data);
    * // likConst = 3.67575413281869;
    * // if m==n
    * //     for k=1:q
    * //         i = indices(k,1);
    * //         j = indices(k,1);
    * //         s11 = S(i,i);
    * //         s12 = S(i,j);
    * //         s22 = S(j,j);
    * //         det = s11*s22-s12*s12;
    * //         x11 = data(i,i);
    * //         x12 = data(i,j);
    * //         x22 = data(j,j);
    * //         ll = 0.5*(likConst + log(det) + (s22*x11 - 2*s12*x12 + s11*x22)/det)/m;
    * //     end
    * // else
    * //     for k=1:q
    * //         i = indices(k,1);
    * //         j = indices(k,1);
    * //         s11 = S(i,i);
    * //         s12 = S(i,j);
    * //         s22 = S(j,j);
    * //         det = s11*s22-s12*s12;
    * //         x11 = data(i)*data(i);
    * //         x12 = data(i)*data(j);
    * //         x22 = data(j)*data(j);
    * //         ll = 0.5*(likConst + log(det) + (s22*x11 - 2*s12*x12 + s11*x22)/det)/m;
    * //     end
    * // end*/
   size_t k, i, j;
   double s11, s12, s22, x11, x12, x22, scale, det;
   double likConst = 3.67575413281869;
   scale = (double)q;
   *ll = 0.0;
   
   if (m==n)
   {
      for (k=0;k<q;k++)
      {
         i = (size_t)indices[k] - 1;
         j = (size_t)indices[k+q] - 1;
         s11 = S[i*m+i];
         s12 = S[i*m+j];
         s22 = S[j*m+j];
         x11 = X[i*m+i];
         x12 = X[i*m+j];
         x22 = X[j*m+j];
         det = s11*s22 - s12*s12;
         *ll += 0.5*(likConst + log(det) + (s22*x11 - 2*s12*x12 + s11*x22)/det)/scale;
      }
   }
   else
   {
      for (k=0;k<q;k++)
      {
         i = (size_t)indices[k] - 1;
         j = (size_t)indices[k+q] - 1;
         s11 = S[i*m+i];
         s12 = S[i*m+j];
         s22 = S[j*m+j];
         x11 = X[i]*X[i];
         x12 = X[i]*X[j];
         x22 = X[j]*X[j];
         det = s11*s22 - s12*s12;
         *ll += 0.5*(likConst + log(det) + (s22*x11 - 2*s12*x12 + s11*x22)/det)/scale;
      }
   }
}

/*
 * dcc_likestream_core.c -
 * Likelihood for estimation of scalar DCC(1,1) multivarate volatility
 * model with GARCH(1,1) conditional variances
 *
 * Based on dcc_likelihood.m but try to make it run faster
 *
 * To be used in the call to fmincon optimizer in dccstreamXX.m
 */

void dcc_likestream_core(double a, double b, const size_t k, const size_t T, double *stdData, double *R, double *backCast, double *llsum) {
   
   /* int count = 0; */
   /* mexPrintf("dcc_likestream_core marker %d \n",count); */
   /* 0 */
   /* count++; */
   double *intercept, aplusb, omega, *Q, *S, *q, *X, *indices;
   double *ll;
   size_t i, j, t, elements, numindices;
   mxArray *mymatrix;
   /* mexPrintf("dcc_likestream_core marker %d \n",count); */
   /* 1 */
   /* count++; */
   
   /*  Compute intercept. */
   aplusb    = a + b;
   omega     = 1 - aplusb;
   if (omega<0)
      omega = 0;
   elements  = k * k;
   intercept = mxCalloc(elements,sizeof(double));
   /* mexPrintf("dcc_likestream_core marker %d \n",count); */
   /* 2 */
   /* count++; */
   for (i=0;i<elements;i++)
   {
      intercept[i] = R[i] * omega;
   }
   /* mexPrintf("dcc_likestream_core marker %d \n",count); */
   /* 3 */
   /* count++; */
   
   /* Specify indices for composite likelihood */
   numindices = k-1;
   indices = mxCalloc((numindices*2),sizeof(double));
   /* mexPrintf("dcc_likestream_core marker %d \n",count); */
   /* 4 */
   /* count++; */
   /* NOTE: TIME UPDATE LOOP ASSUMES THIS PARTICULAR PATTERN OF INDICES */
   /* IF INDICES ARE MODIFIED, THEN TIME UPDATE LOOP MUST BE REWRITTEN */
   /* THIS SAVES TIME WITH COMPOSITE LIKELIHOOD */
   for (i=0;i<numindices;i++)
   {
      indices[i]            = (double) i + 1;
      indices[i+numindices] = (double) i + 2;
   }
   /* mexPrintf("dcc_likestream_core marker %d \n",count); */
   /* 5 */
   /* count++; */
   
   /* create arrays */
   /* mexPrintf("dcc_likestream_core marker %d \n",count); */
   /* 6 */
   /* count++; */
   Q = mxCalloc(elements,sizeof(double));
   S = mxCalloc(elements,sizeof(double));
   q = mxCalloc(k,sizeof(double));
   X = mxCalloc(elements,sizeof(double));
   /* mexPrintf("dcc_likestream_core marker %d \n",count); */
   /* 7 */
   /* count++; */
   
   /* initialize likelihoods */
   /* mexPrintf("dcc_likestream_core marker %d \n",count); */
   /* 8 */
   /* count++; */
   *llsum  = 0.0;
   mymatrix = mxCreateDoubleMatrix(1,1,mxREAL);
   ll = mxGetPr(mymatrix);
   /* mexPrintf("dcc_likestream_core marker %d \n",count); */
   /* 9 */
   /* count++; */
   
   /* loop over dates */
   /* mexPrintf("%s \n","I am starting the loop now and so far so good."); */
   for (t=0;t<T;t++)
   {
      
      /* mexPrintf("%d \n",t); */
      /* update conditional correlation matrix */
      /* mexPrintf("%s \n","update conditional correlation matrix"); */
      /* NOTE: TIME UPDATE LOOP ASSUMES A SPECIFIC PATTERN OF INDICES */
      /* IF YOU MODIFY INDICES, YOU MUST ADJUST THIS LOOP ACCORDINGLY */
      /* THIS IS TO SAVE TIME WITH COMPOSITE LIKELIHOOD */
      if (t==0)
      {
         Q[0] = intercept[0] + (aplusb*backCast[0]);
         for (i=1;i<k;i++) 
         {
            for (j=i-1;j<=i;j++)
            {
               Q[i+(k*j)] = intercept[i+(k*j)] + (aplusb*backCast[i+(k*j)]);
            }
         }
      }
      else
      {
         Q[0]=intercept[0] + (a*stdData[(elements*(t-1))]) + (b*Q[0]);
         for (i=1;i<k;i++) 
         {
            for (j=i-1;j<=i;j++)
            {
               Q[i+(k*j)]=intercept[i+(k*j)] + (a*stdData[i+(k*j)+(elements*(t-1))]) + (b*Q[i+(k*j)]);
            }
         }
      }
      
      /* extract the square root of the diagonal elements */
      /* mexPrintf("%s \n","extract the square root of the diagonal elements"); */
      for (i=0;i<k;i++)
      {
         q[i] = sqrt(Q[i+(k*i)]);
      }
      
      /* set the diagonal of the conditional correlation matrix to one */
      /* mexPrintf("%s \n","set the diagonal of the conditional correlation matrix to one"); */
      S[0] = 1;
      for (i=1;i<k;i++)
      {
         /* mexPrintf("i = %d \n",i); */
         for (j=i-1;j<=i;j++)
         {
            /* mexPrintf("i = %d , ",i);              */
            /* mexPrintf("j = %d , ",j);              */
            /* mexPrintf("k = %d , ",k);              */
            /* mexPrintf("i+(k*j) = %d \n ",i+(k*j)); */
            S[i+(k*j)] = Q[i+(k*j)]/(q[i]*q[j]);
         }
      }
      
      /* extract current day of standardized residuals */
      /* mexPrintf("%s \n","extract current day of standardized residuals"); */
      X[0] = stdData[elements*t];
      for (i=1;i<k;i++)
      {
         for (j=i-1;j<=i;j++)
         {
            X[i+(k*j)] = stdData[i+(k*j)+(elements*t)];
         }
      }
      
      /* call Kevin Sheppard's function to compute composite likelihood */
      /* mexPrintf("%s \n","call Kevin Sheppards function to compute composite likelihood"); */
      composite_likelihood_core(S, X, indices, numindices, k, k, ll);
      
      /* sum up the individual likelihoods */
      /* mexPrintf("%s \n","sum up the individual likelihoods"); */
      *llsum += *ll;
   }
   /* mexPrintf("%s \n","I am finishing the loop now."); */
   /* mexPrintf("dcc_likestream_core marker %d \n",count); */
   /* 10 */
   /* count++; */
   
   /* clear up the memory */
   mxDestroyArray(mymatrix);
   mxFree(intercept);
   mxFree(indices);
   mxFree(Q);
   mxFree(S);
   mxFree(q);
   mxFree(X);
   /* mexPrintf("dcc_likestream_core marker %d \n",count); */
   /* 11 */
   /* count++; */
}

/* The gateway routine */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
   
   /* Note: first input (parameters) must be a VERTICAL 2*1 vector */
   
   double *parameters, *stdData, *R, *backCast, a, b;
   int j;
   size_t k, T;
   const mwSize *dim_array;
   double *llsum;
   clock_t t;
   double time_taken;

   t = clock();
   
   /*  Check for proper number of arguments. */
   if (nrhs!=4)
      mexErrMsgTxt("Three inputs required.");
   if (nlhs!=1)
      mexErrMsgTxt("One output required.");
   
   /* Check data type of input arguments */
   for (j=0; j<4; j++)
   {
      if (!(mxIsDouble(prhs[j])) || (mxIsComplex(prhs[j]))) {
         mexErrMsgTxt("All input arrays must be of type real double.");
      }
   }
   
   /* Check inputs have the required numbers of dimensions */
   if (mxGetNumberOfDimensions(prhs[0])!=2)
      mexErrMsgTxt("First input should be 2-dimensional matrix");
   if (mxGetNumberOfDimensions(prhs[1])!=3)
      mexErrMsgTxt("Second input should be 3-dimensional array");
   if (mxGetNumberOfDimensions(prhs[2])!=2)
      mexErrMsgTxt("Third input should be 2-dimensional matrix");
   if (mxGetNumberOfDimensions(prhs[3])!=2)
      mexErrMsgTxt("Fourth input should be 2-dimensional matrix");
   
   /* Get the number of variables and the sample size */
   k = (size_t)mxGetM(prhs[2]);
   dim_array = mxGetDimensions(prhs[1]);
   T = (size_t)dim_array[2];
   /* mexPrintf("0 %d \n",dim_array[0]); */
   /* mexPrintf("1 %d \n",dim_array[1]); */
   /* mexPrintf("2 %d \n",dim_array[2]); */
   /* mexPrintf("%s %d \n","T=",T); */
   /* mexPrintf("%s %d \n","k=",k); */
   
   /* Check the dimensions of the inputs are correct*/
   if (mxGetM(prhs[0])!=1)
      mexErrMsgTxt("First input should be horizontal vector.");
   if (mxGetN(prhs[0])!=2)
      mexErrMsgTxt("First input should vector with two elements.");
   if (mxGetN(prhs[2])!=k)
      mexErrMsgTxt("Third input should be a square matrix.");
   if (mxGetM(prhs[3])!=k)
      mexErrMsgTxt("4th input should have same number of rows as 3rd.");
   if (mxGetN(prhs[3])!=k)
      mexErrMsgTxt("4th input must have same number of columns as 3rd");
   
   /* Check dimensions of the stdData matrix */
   if (dim_array[0]!=k)
      mexErrMsgTxt("2nd input should have same number of rows as 3rd.");
   if (dim_array[1]!=k)
      mexErrMsgTxt("2nd input must have same number of columns as 3rd");
   
   /*  Create a pointer to the input matrices . */
   parameters = mxGetPr(prhs[0]);
   stdData    = mxGetPr(prhs[1]);
   R          = mxGetPr(prhs[2]);
   backCast   = mxGetPr(prhs[3]);
   
   /* Extract the ARCH/GARCH parameters */
   a = parameters[0];
   b = parameters[1];
   
   /*  Set the output pointer to the output matrix. */
   plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
   
   /*  Create a C pointer to a copy of the output matrix. */
   llsum = mxGetPr(plhs[0]);
   
   /*  Call the C subroutine. */
   /* mexPrintf("%s \n","Calling the C subroutine."); */
   dcc_likestream_core(a,b,k,T,stdData,R,backCast,llsum);
   
   /* sanity checks */
   if (!mxIsFinite(*llsum))
      *llsum=10000000.0;
   
   /* print inputs and outputs */
   t = clock() - t;
   time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
   if (1==2) {
      mexPrintf("a = %16.14f, ", a);
      mexPrintf("b = %16.14f, ", b);
      mexPrintf("likelihood = %16.14f\n", *llsum);
      printf("dcc_likestream05 in C took %f seconds to execute \n", time_taken);
   }   
}
