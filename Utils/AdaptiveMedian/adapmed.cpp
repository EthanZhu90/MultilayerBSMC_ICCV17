#include <vector>
#include <algorithm>
#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[])
{
  // Check for proper number of arguments
  if (nrhs != 2)
    mexErrMsgIdAndTxt("Mex:adapmed:nrhs", "Only two input expected.");

  if (nlhs != 1)
    mexErrMsgIdAndTxt("Mex:ldoftrack:nlhs", "One output required.");
  
  if (mxGetNumberOfElements(prhs[1]) != mxGetN(prhs[0]))
    mexErrMsgIdAndTxt("Mex:adapmed", "Number of counts should match number of columns.");

  mwSignedIndex m,n;
  m = mxGetM(prhs[0]);
  const mwSize* dims = mxGetDimensions(prhs[0]);
  n = dims[1];

  mxArray* arr = mxDuplicateArray(prhs[0]);
  double* data = (double*) mxGetData(arr);
  double* counts = (double*) mxGetData(prhs[1]);

  // Create output array
  plhs[0] = mxCreateDoubleMatrix(n, 1,mxREAL);
  double* out = (double*)mxGetData(plhs[0]);

  double* p = data;
  double* c = counts;

  for(mwSignedIndex i=0; i<n; ++i) {
    // Sort elements in input array
    int cc = (int)*c;
    if (cc > m) {
      mexErrMsgIdAndTxt("Mex:adapmed", "Out of bounds count!");
    }

    if (cc > 2) {
      std::sort(p, p+cc);
    }

    // Extract the median
    if (cc == 0) {
      *out = 0;
    } else if (cc % 2 == 0) {
      *out = (p[cc/2 - 1] + p[cc/2]) / 2;
    } else {
      *out = p[cc/2];
    }

    p+= m;
    c++;
    out++;
  }
}
