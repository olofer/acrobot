/*
 * Octave/C++ ODE-type function for state variable z = (theta1, theta2, theta1dot, theta2dot).
 * Parameter struct P must contain {M1, M2, I1, I2, L1, L2, g, muA, muB}
 *
 * USAGE: dz = acrobot_odefun(t, z, P, u);            % u optional
 * BUILD: see "acrobot_rebuild_test.m"
 *
 * ALT. USAGE: [A, B] = acrobot_odefun(ep, z, P, u);   % u optional, return linearization @ (z, u)
 *
 */

#include "mex.h"
#include <cstring>
#include <cmath>

#include "../acrobot.hpp"
#include "acrobot_mex_utils.hpp"

void mexFunction(int nlhs, 
                 mxArray** plhs, 
                 int nrhs, 
                 const mxArray** prhs)
{
  const int arg_index_time = 0;
  const int arg_index_state = 1;
  const int arg_index_params = 2;
  const int arg_index_torque = 3;

  if ((nrhs < 3 || nrhs > 4) || nlhs > 2) {
    mexErrMsgTxt("USAGE: dz = acrobot_odefun(t, z, P [, u]);");
  }

  if (!isDoubleRealScalar(prhs[arg_index_time])) {
    mexErrMsgTxt("Expects t to be real scalar");
  }

  if (!(isDoubleRealVector(prhs[arg_index_state]) && 
        mxGetNumberOfElements(prhs[arg_index_state]) == 4))
  {
    mexErrMsgTxt("Expects z to be a real vector with 4 elements");
  }

  if (!mxIsStruct(prhs[arg_index_params])) {
    mexErrMsgTxt("Expects P to be a (scalar) struct");
  }

  const double time = getDoubleScalar(prhs[arg_index_time]);
  const double* ptrz = (const double*) mxGetPr(prhs[arg_index_state]); 

  acrobot::params P;
  set_from_struct(P, prhs[arg_index_params]);
  if (!P.sanityChecksOut()) {
    mexErrMsgTxt("P does not hold a complete and valid set of parameters");
  }

  if (nrhs == 4) {
    if (!isDoubleRealScalar(prhs[arg_index_torque])) {
      mexErrMsgTxt("Expects u to be real scalar");
    }
    const double torque = getDoubleScalar(prhs[arg_index_torque]);
    P.u = torque;
  }

  if (nlhs == 1) {
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[arg_index_state]), mxGetN(prhs[arg_index_state]), mxREAL);
    double* ptrdz = (double *) mxGetPr(plhs[0]);

    acrobot::calculate_dotted_state(ptrdz, time, ptrz, &P);
    return;
  }

  if (nlhs != 2) {
    mexErrMsgTxt("Must have exactly 1 or 2 outputs requested");
  }

  // NOTE: use time as finite difference epsilon (it is never actually used anyway)
  const double ep = time;
  if (ep <= 0.0 || ep > 0.5) {
    mexErrMsgTxt("FD epsilon is out of (sane) range (it is specified by the t parameter)");
  }

  plhs[0] = mxCreateDoubleMatrix(4, 4, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(4, 1, mxREAL);

  acrobot::linearize_fd((double*) mxGetPr(plhs[0]),
                        (double*) mxGetPr(plhs[1]),
                        ptrz,
                        P.u,
                        &P,
                        ep);
  return;
}
