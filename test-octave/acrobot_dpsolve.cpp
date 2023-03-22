/*
 * USAGE: [V, A] = acrobot_dpsolve(P, itrs);
 *
 */

// TODO: need the develop the neccessary equation/stencil for the evolution rule to be used
//       parameters include "dt" and "sigma" and the dynamics "dotted state" function
// TODO: also figure out correct treatment of periodic state variables and absorbing state variables

#include "mex.h"
#include <cstring>
#include <cmath>

#include "../acrobot.hpp"
#include "acrobot_mex_utils.hpp"

#include <vector>

const double _one_pi = 3.14159265358979323846;

class acrobotDP
{
public:
  acrobotDP(int nth1, 
            int nth2, 
            int ndot1, 
            int ndot2,
            double maxdot1 = 1.0,
            double maxdot2 = 1.0) : 
  grid_th1(nth1), 
  grid_th2(nth2), 
  grid_th1d(ndot1), 
  grid_th2d(ndot2),
  value(nth1 * nth2 * ndot1 * ndot2), 
  action(nth1 * nth2 * ndot1 * ndot2), 
  dims{nth1, nth2, ndot1, ndot2}
  {
    for (int i = 0; i < nth1; i++) {
      grid_th1[i] = (2.0 * _one_pi * i) / nth1; 
    }
    for (int i = 0; i < nth2; i++) {
      grid_th2[i] = (2.0 * _one_pi * i) / nth2; 
    }
    for (int i = 0; i < ndot1; i++) {
      grid_th1d[i] = -maxdot1 + (2.0 * maxdot1 * i) / (ndot1 - 1);
    }
    for (int i = 0; i < ndot2; i++) {
      grid_th2d[i] = -maxdot2 + (2.0 * maxdot2 * i) / (ndot2 - 1);
    }
  }

  ~acrobotDP() {}

  int getN1() const { return grid_th1.size(); }
  int getN2() const { return grid_th2.size(); }
  int getN1D() const { return grid_th1d.size(); }
  int getN2D() const { return grid_th2d.size(); }
  int size() const { return value.size(); }
  int dim(int i) const { return dims[i]; }
  int sub2ind(const int* i) const { return sub2ind(i[0], i[1], i[2], i[3]); }
  int sub2ind(int i0, int i1, int i2, int i3) const {
    //return i0 + dims[0] * i1 + dims[0] * dims[1] * i2 + dims[0] * dims[1] * dims[2] * i3;
    return i0 + dims[0] * (i1 + dims[1] * (i2 + dims[2] * i3));
  }
  void ind2sub(int idx, int* i) const {
    i[0] = idx % dims[0];
    idx = (idx - i[0]) / dims[0];
    i[1] = idx % dims[1];
    idx = (idx - i[1]) / dims[1];
    i[2] = idx % dims[2];
    idx = (idx - i[2]) / dims[2];
    i[3] = idx; 
  }
  bool check_ind2sub2ind() const {
    int indices[4];
    for (int k = 0; k < size(); k++) {
      ind2sub(k, indices);
      if (k != sub2ind(indices))
        return false;
    }
    return true;
  }

  double errsq(const acrobot::params* P, // given P->u provided control value
               const double* xsource, 
               const double* xtarget, 
               double dt,
               const double* sigma) const
  {
    // Use trapezoidal transcription equation to assess the "fitness" of the source and target states w.r.t. dynamics
    double deltx[4] = {xtarget[0] - xsource[0], xtarget[1] - xsource[1], xtarget[2] - xsource[2], xtarget[3] - xsource[3]};
    if (deltx[0] < -_one_pi) deltx[0] += 2.0 * _one_pi;
    if (deltx[0] > _one_pi)  deltx[0] -= 2.0 * _one_pi;
    if (deltx[1] < -_one_pi) deltx[1] += 2.0 * _one_pi;
    if (deltx[1] > _one_pi)  deltx[1] -= 2.0 * _one_pi;

    const double dummytime = 0.0;
    double dotx0[4] = {0.0, 0.0, 0.0, 0.0};
    acrobot::calculate_dotted_state(dotx0, dummytime, xsource, P);
    double dotx1[4] = {0.0, 0.0, 0.0, 0.0};
    acrobot::calculate_dotted_state(dotx1, dummytime, xtarget, P);

    double errsq = 0.0;
    for (int i = 0; i < 4; i++) {
      const double wi = deltx[i] - 0.5 * dt * (dotx0[i] + dotx1[i]);
      errsq += wi * wi / (sigma[i] * sigma[i] * dt);
    }
    return errsq;
  }

  // FIXME: member function that evaluates all control options across all neighbor target lattice points from a given source lattice point

  int initilize() {
    // set terminal cost ; the value is just the number of time-steps to-go until done
    // the action array A can signal whether there has been any value assigned at a cell
    return 0;
  }

  int update() {
    // run a single scan across the cells; update values where it is possible to do so; return number of values edited
    return 0;
  }

  // FIXME: missing member functions to "apply" the control law: i.e. control = evaluate(state).

private:
  std::vector<double> grid_th1;
  std::vector<double> grid_th2;
  std::vector<double> grid_th1d;
  std::vector<double> grid_th2d;

  int dims[4];

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

  acrobotDP acbdp(36, 34, 33, 35);

  mexPrintf("num. cells = %i\n", acbdp.size());
  mexPrintf("dims = %i,%i,%i,%i\n", acbdp.dim(0), acbdp.dim(1), acbdp.dim(2), acbdp.dim(3));

  if (!acbdp.check_ind2sub2ind()) {
    mexErrMsgTxt("Self-check of ind2sub, sub2ind failed");
  }
  
  acbdp.initilize();

  const mwSize dims[4] = {acbdp.getN1(), acbdp.getN2(), acbdp.getN1D(), acbdp.getN2D()};
  plhs[0] = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);

  double* V = (double*) mxGetPr(plhs[0]);
  double* A = (double*) mxGetPr(plhs[1]);

  for (int k = 0; k < acbdp.size(); k++) {
    int indices[4];
    acbdp.ind2sub(k, indices);
    V[acbdp.sub2ind(indices[0], indices[1], indices[2], indices[3])] = (double) k;
    A[acbdp.sub2ind(indices[0], indices[1], indices[2], indices[3])] = (double) (acbdp.size() - k - 1);
  }

  return;
}
