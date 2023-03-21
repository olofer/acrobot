/*
 * USAGE: [V, A] = acrobot_dpsolve(P, itrs);
 *
 */

#include "mex.h"
#include <cstring>
#include <cmath>

#include "../acrobot.hpp"
#include "acrobot_mex_utils.hpp"

#include <vector>

class acrobotDP
{
public:
  acrobotDP(int nth1, 
            int nth2, 
            int ndot1, 
            int ndot2) : 
  grid_th1(nth1), grid_th2(nth2), grid_th1d(ndot1), grid_th2d(ndot2),
  value(nth1 * nth2 * ndot1 * ndot2), action(nth1 * nth2 * ndot1 * ndot2)
  {
    // fill in the grid locations
    // allocate the cell 4D space; and initialize it to zeros or "nans"
    // prepare convenience indexing factors
  }

  ~acrobotDP() {}

  int getN1() const { return grid_th1.size(); }
  int getN2() const { return grid_th2.size(); }
  int getN1D() const { return grid_th1d.size(); }
  int getN2D() const { return grid_th2d.size(); }

  int size() const { return value.size(); }

  int initilize() {
    // set terminal cost ; the value is just the number of time-steps to-go until done
    return 0;
  }

  int update() {
    // run a single scan across the cells; update values where it is possible to do so; return number of values edited
    return 0;
  }

private:
  std::vector<double> grid_th1;
  std::vector<double> grid_th2;
  std::vector<double> grid_th1d;
  std::vector<double> grid_th2d;

  std::vector<float> value;
  std::vector<int8_t> action;
};

void mexFunction(int nlhs, 
                 mxArray** plhs, 
                 int nrhs, 
                 const mxArray** prhs)
{
  const int arg_index_params = 0;
  const int arg_index_iterations = 1;

  if (nrhs != 2 || nlhs != 2) {
    mexErrMsgTxt("USAGE: [V, A] = acrobot_dpsolve(P, itrs);");
  }

  if (!isDoubleRealScalar(prhs[arg_index_iterations])) {
    mexErrMsgTxt("Expects itrs to be real scalar");
  }

  const int itrs = (int) std::round(getDoubleScalar(prhs[arg_index_iterations]));

  if (itrs <= 0) {
    mexErrMsgTxt("itrs >= 1 required");
  }

  if (!mxIsStruct(prhs[arg_index_params])) {
    mexErrMsgTxt("Expects P to be a (scalar) struct");
  }

  acrobot::params P;
  set_from_struct(P, prhs[arg_index_params]);
  if (!P.sanityChecksOut()) {
    mexErrMsgTxt("P does not hold a complete and valid set of parameters");
  }

  //plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[arg_index_state]), mxGetN(prhs[arg_index_state]), mxREAL);
  //double* ptrdz = (double *) mxGetPr(plhs[0]);
  //acrobot::calculate_dotted_state(ptrdz, time, ptrz, &P);

  acrobotDP acbdp(9, 10, 11, 12);

  mexPrintf("num. cells = %i\n", acbdp.size());

  const mwSize dims[4] = {acbdp.getN1(), acbdp.getN2(), acbdp.getN1D(), acbdp.getN2D()};
  plhs[0] = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);

  return;
}
