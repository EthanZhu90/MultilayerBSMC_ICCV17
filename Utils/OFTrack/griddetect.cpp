#include <vector>
#include <algorithm>
#include <ctime>
#include "CTensor.h"
#include "CFilter.h"
#include "mex.h"

class CTrack {
public:
  CTrack() {}
  float mx,my;             // current position of the track
  float mox,moy;                          // original starting point of the track
};

int mStep;
int mXSize,mYSize;

// computeCorners --------------------------------------------------------------
void computeCorners(CTensor<float>& aImage, CMatrix<float>& aCorners, float aRho) {
  aCorners.setSize(aImage.xSize(),aImage.ySize());
  int aXSize = aImage.xSize();
  int aYSize = aImage.ySize();
  int aSize = aXSize*aYSize;
  // Compute gradient
  CTensor<float> dx(aXSize,aYSize,aImage.zSize());
  CTensor<float> dy(aXSize,aYSize,aImage.zSize());
  CDerivative<float> aDerivative(3);
  NFilter::filter(aImage,dx,aDerivative,1,1);
  NFilter::filter(aImage,dy,1,aDerivative,1);
  // Compute second moment matrix
  CMatrix<float> dxx(aXSize,aYSize,0);
  CMatrix<float> dyy(aXSize,aYSize,0);
  CMatrix<float> dxy(aXSize,aYSize,0);
  int i2 = 0;
  for (int k = 0; k < aImage.zSize(); k++)
    for (int i = 0; i < aSize; i++,i2++) {
      dxx.data()[i] += dx.data()[i2]*dx.data()[i2];
      dyy.data()[i] += dy.data()[i2]*dy.data()[i2];
      dxy.data()[i] += dx.data()[i2]*dy.data()[i2];
    }
  // Smooth second moment matrix
  NFilter::recursiveSmoothX(dxx,aRho);
  NFilter::recursiveSmoothY(dxx,aRho);
  NFilter::recursiveSmoothX(dyy,aRho);
  NFilter::recursiveSmoothY(dyy,aRho);
  NFilter::recursiveSmoothX(dxy,aRho);
  NFilter::recursiveSmoothY(dxy,aRho);
  // Compute smallest eigenvalue
  for (int i = 0; i < aSize; i++) {
    float a = dxx.data()[i];
    float b = dxy.data()[i];
    float c = dyy.data()[i];
    float temp = 0.5*(a+c);
    float temp2 = temp*temp+b*b-a*c;
    if (temp2 < 0.0f) aCorners.data()[i] = 0.0f;
    else aCorners.data()[i] = temp-sqrt(temp2);
  }
}

// Code fragment from Pedro Felzenszwalb  ---------------
// http://people.cs.uchicago.edu/~pff/dt/
void dt(CVector<float>& f, CVector<float>& d, int n) {
  d.setSize(n);
  int *v = new int[n];
  float *z = new float[n+1];
  int k = 0;
  v[0] = 0;
  z[0] = -10e20;
  z[1] = 10e20;
  for (int q = 1; q <= n-1; q++) {
    float s  = ((f[q]+q*q)-(f(v[k])+v[k]*v[k]))/(2*(q-v[k]));
    while (s <= z[k]) {
      k--;
      s  = ((f[q]+q*q)-(f(v[k])+v[k]*v[k]))/(2*(q-v[k]));
    }
    k++;
    v[k] = q;
    z[k] = s;
    z[k+1] = 10e20;
  }
  k = 0;
  for (int q = 0; q <= n-1; q++) {
    while (z[k+1] < q)
      k++;
    int help = q-v[k];
    d(q) = help*help + f(v[k]);
  }
  delete[] v;
  delete[] z;
}
// ------------------------------------------------------

// euclideanDistanceTransform
void euclideanDistanceTransform(CMatrix<float>& aMatrix) {
  int aXSize = aMatrix.xSize();
  int aYSize = aMatrix.ySize();
  CVector<float> f(NMath::max(aXSize,aYSize));
  // Transform along columns
  for (int x = 0; x < aXSize; x++) {
    for (int y = 0; y < aYSize; y++)
      f(y) = aMatrix(x,y);
    CVector<float> d;
    dt(f,d,aYSize);
    for (int y = 0; y < aYSize; y++)
      aMatrix(x,y) = d(y);
  }
  // Transform along rows
  for (int y = 0; y < aYSize; y++) {
    int aOffset = y*aXSize;
    for (int x = 0; x < aXSize; x++)
      f(x) = aMatrix.data()[x+aOffset];
    CVector<float> d;
    dt(f,d,aXSize);
    for (int x = 0; x < aXSize; x++)
      aMatrix.data()[x+aOffset] = d(x);
  }
}

void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[])
{
  std::vector<CTrack> mTracks;
  std::vector<CTrack> mTracks2;

  // Check for proper number of arguments
  if (nrhs > 3)
    mexErrMsgIdAndTxt("Mex:ldoftrack:nrhs", "At most three inputs.");

  if (nlhs != 1)
    mexErrMsgIdAndTxt("Mex:ldoftrack:nlhs", "One output required.");
  
  // Input 1 is the first image : it should be a matrix of unsinged int 8 
  if (!mxIsUint8(prhs[0]))
    mexErrMsgTxt("im must be a full matrix of uint8 values.");

  if (nrhs == 3)
    if (mxGetM(prhs[2]) != 2)
      mexErrMsgTxt("x must be 2xN.");

  mStep = (int) mxGetScalar(prhs[1]);

  // We have the number of rows and columns reversed since matlab uses column major 
  mwSignedIndex m,n;
  mwSignedIndex numFeatures = 0;
  m = mxGetM(prhs[0]);
  const mwSize* dims = mxGetDimensions(prhs[0]);
  n = dims[1];
  
  unsigned char* im1 = (unsigned char*) mxGetData(prhs[0]);
  double* x = NULL;
  if (nrhs == 3) {
    numFeatures = mxGetN(prhs[2]);   // existing point
    x = (double*)mxGetData(prhs[2]);
  }

  srand(time(NULL));

  // Copy images into row-major format.
  CTensor<float> aImage1(n, m, 3);
  
  for (int ac=0; ac<3; ++ac) {
    for (int ax=0; ax< n; ++ax) {
      for (int ay=0; ay<m; ++ay) {
	aImage1(ax,ay,ac) = *(im1 + ac * m * n + ax * m + ay);
      }
    }
  }

  mTracks.clear();
  // Add existing tracks
  for (int i=0;i<numFeatures;++i) {
    mTracks.push_back(CTrack());
    CTrack& track = mTracks.back();
    track.mox = x[2*i] - 1.0;
    track.moy = x[2*i + 1] - 1.0;
  }

  mXSize = aImage1.xSize();
  mYSize = aImage1.ySize();
  CMatrix<float> aCorners;
  CMatrix<float> aCovered(mXSize,mYSize);
  int aSize = mXSize*mYSize;
  // Smooth first image (can be removed from optical flow computation then)
  NFilter::recursiveSmoothX(aImage1,0.8f);
  NFilter::recursiveSmoothY(aImage1,0.8f);

  // Tracking
  // Mark areas sufficiently covered by tracks
  aCovered = 1e20;
  if (mTracks.size() > 0) {
    for (unsigned int i = 0; i < mTracks.size(); i++) {
      if (mTracks[i].mox < 0 || mTracks[i].mox > (mXSize-1) || mTracks[i].moy < 0 || mTracks[i].moy > (mYSize-1)) {
	//mexPrintf("Out of bounds points detected in input!\n");
	continue;
      }
      aCovered((int)mTracks[i].mox,(int)mTracks[i].moy) = 0.0f;
    }
    euclideanDistanceTransform(aCovered);
  }

  // Set up new tracking points in uncovered areas
  computeCorners(aImage1,aCorners,3.0f);
  float aCornerAvg = aCorners.avg();
  
  for (int ay = 4; ay < mYSize-4; ay+=mStep) {
    for (int ax = 4; ax < mXSize-4; ax+=mStep) {
      if (aCovered(ax,ay) < mStep*mStep) continue;
      float distToImageBnd = exp(-0.1*NMath::min(NMath::min(NMath::min(ax,ay),mXSize-ax),mYSize-ay));
     // if (aCorners(ax,ay) < 1.0* (aCornerAvg*(0.1+distToImageBnd))) continue;  // orginal 
      
      float thresholdAvg = 1.0* (aCornerAvg*(0.1+distToImageBnd));  // revised by Ethan 11/16
      if (aCorners(ax,ay) < thresholdAvg)
      {
          float prob = 0.5*(cos((aCorners(ax,ay)/thresholdAvg - 1) * M_PI ) + 1);  
          if(rand()/(RAND_MAX + 1.0) > prob) continue; 
      }



      if (aCorners(ax,ay) < 1.0*(1.0f+distToImageBnd)) continue;
      mTracks2.push_back(CTrack());
      CTrack& newTrack = mTracks2.back();
      newTrack.mox = ax;
      newTrack.moy = ay;
    }
  }
 
  // Create output array
  plhs[0] = mxCreateDoubleMatrix(2, mTracks2.size(),mxREAL);
  double* xout = (double*)mxGetData(plhs[0]);

  // Copy results back  to output array
  for(int i=0;i<mTracks2.size();i++) {
    xout[2*i] = mTracks2[i].mox + 1;
    xout[2*i+1] = mTracks2[i].moy + 1;
  }

}


