#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include "CTensor.h"
#include "CFilter.h"
#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[])
{
  // Check for proper number of arguments
  if (nrhs !=3)
	  mexErrMsgIdAndTxt("Mex:ldoftrack:nrhs", "Five inputs required.");

  if (nlhs != 2)
	  mexErrMsgIdAndTxt("Mex:ldoftrack:nlhs", "At least one output required.");
  
  // Input 1 is the first image : it should be a matrix of unsinged int 8 
  if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]))
	  mexErrMsgTxt("flow1 and flow2 must be a full matrix of double values.");
  if (mxGetM(prhs[0]) != mxGetM(prhs[1]) || mxGetN(prhs[0]) != mxGetN(prhs[1]))
	  mexErrMsgTxt("flow1 and flow2 must have the same dimensions.");

  if (mxGetM(prhs[2]) != 2)
	  mexErrMsgTxt("x must be 2xN.");

  // We have the number of rows and columns reversed since matlab uses column major 
  mwSignedIndex m,n, numFeatures;
  m = mxGetM(prhs[0]);
  const mwSize* dims = mxGetDimensions(prhs[0]);
  n = dims[1];
  numFeatures = mxGetN(prhs[2]);

  double* forward = (double*) mxGetData(prhs[0]);
  double* backward = (double*) mxGetData(prhs[1]);
  double* x = (double*)mxGetData(prhs[2]);

   // Create output array
  plhs[0] = mxCreateDoubleMatrix(2, numFeatures,mxREAL);
  plhs[1] = mxCreateLogicalMatrix(1, numFeatures);
  double* xout = (double*)mxGetData(plhs[0]);
  mxLogical* tracked = (mxLogical*) mxGetData(plhs[1]);
 
  // Copy flow to tensor format
  CTensor<float> aForward(n,m,2);
  CTensor<float> aBackward(n,m,2);
  for (int ac=0; ac<2; ++ac) {
    for (int ax=0; ax< n; ++ax) {
      for (int ay=0; ay<m; ++ay) {
	aForward(ax,ay,ac) = *(forward + ac * m * n + ax * m + ay);
	aBackward(ax,ay,ac) = *(backward + ac * m * n + ax * m + ay);
      }
    }
  }

  int mXSize = aForward.xSize();
  int mYSize = aForward.ySize();
  int aSize = mXSize*mYSize;

  // Check consistency of forward flow via backward flow
  CMatrix<float> aUnreliable(mXSize,mYSize,0);
  CTensor<float> dx(mXSize,mYSize,2);
  CTensor<float> dy(mXSize,mYSize,2);
  CDerivative<float> aDev(3);
  NFilter::filter(aForward,dx,aDev,1,1);
  NFilter::filter(aForward,dy,1,aDev,1);
  CMatrix<float> aMotionEdge(mXSize,mYSize,0);
  for (int i = 0; i < aSize; i++) {
    aMotionEdge.data()[i] += dx.data()[i]*dx.data()[i];
    aMotionEdge.data()[i] += dx.data()[aSize+i]*dx.data()[aSize+i];
    aMotionEdge.data()[i] += dy.data()[i]*dy.data()[i];
    aMotionEdge.data()[i] += dy.data()[aSize+i]*dy.data()[aSize+i];
  }
  for (int ay = 0; ay < aForward.ySize(); ay++) {
    for (int ax = 0; ax < aForward.xSize(); ax++) {
      float bx = ax+aForward(ax,ay,0);
      float by = ay+aForward(ax,ay,1);
      int x1 = (int)bx;
      int y1 = (int)by;
      int x2 = x1+1;
      int y2 = y1+1;
      if (x1 < 0 || x2 >= mXSize || y1 < 0 || y2 >= mYSize) { aUnreliable(ax,ay) = 1.0f; continue;}
      float alphaX = bx-x1; float alphaY = by-y1;
      float a = (1.0-alphaX)*aBackward(x1,y1,0)+alphaX*aBackward(x2,y1,0);
      float b = (1.0-alphaX)*aBackward(x1,y2,0)+alphaX*aBackward(x2,y2,0);
      float u = (1.0-alphaY)*a+alphaY*b;
      a = (1.0-alphaX)*aBackward(x1,y1,1)+alphaX*aBackward(x2,y1,1);
      b = (1.0-alphaX)*aBackward(x1,y2,1)+alphaX*aBackward(x2,y2,1);
      float v = (1.0-alphaY)*a+alphaY*b;
      float cx = bx+u;
      float cy = by+v;
      float u2 = aForward(ax,ay,0);
      float v2 = aForward(ax,ay,1);
      if (((cx-ax)*(cx-ax)+(cy-ay)*(cy-ay)) >= 0.01*(u2*u2+v2*v2+u*u+v*v)+0.5f) { aUnreliable(ax,ay) = 1.0f; continue;}
      if (aMotionEdge(ax,ay) > 0.01*(u2*u2+v2*v2)+0.002f) { aUnreliable(ax,ay) = 1.0f; continue;}
    }
  }

  for (unsigned int i = 0; i < numFeatures; i++) {
    float ax,ay,oldvar;
    ax = x[2*i] - 1; ay = x[2*i+1] - 1;
    int iax = lroundf(ax);
    int iay = lroundf(ay);

    // Disregard unreliable points in bidirectional flow field
    if (aUnreliable(iax,iay) > 0) tracked[i] = false;
    else {
      float bx = ax+aForward(iax,iay,0);
      float by = ay+aForward(iax,iay,1);
      int ibx = lroundf(bx);
      int iby = lroundf(by);

      // Disregard points moving out of the frame
      if (ibx < 0 || iby < 0 || ibx >= mXSize || iby >= mYSize) tracked[i] = false;
      else {
	tracked[i] = true;
	xout[2*i] = bx + 1;
	xout[2*i+1] = by + 1;
      }
    }
  }
}


