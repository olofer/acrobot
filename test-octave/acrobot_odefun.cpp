/*
 * Octave/C++ ODE-type function for state variable z = (theta1, theta2, theta1dot, theta2dot).
 * Parameter struct P must contain {M1, M2, I1, I2, L1, L2, g, muA, muB}
 *
 * USAGE: dz = acrobot_odefun(t, z, P, u);            % u optional
 * BUILD: see "acrobot_rebuild_test.m"
 *
 */

#include "mex.h"
#include <cstring>
#include <cmath>
#include "../acrobot.hpp"

double getDoubleScalar(const mxArray* a);
bool isDoubleRealScalar(const mxArray* a);
bool isDoubleRealVector(const mxArray* a);
bool isDoubleRealMatrix(const mxArray* a);

void set_from_struct(acrobot::params& P, 
                     const mxArray* a)
{
  mxArray* af = nullptr;
  af = mxGetField(a, 0, "M1"); if (isDoubleRealScalar(af))  P.m1 = getDoubleScalar(af);
  af = mxGetField(a, 0, "I1"); if (isDoubleRealScalar(af))  P.I1 = getDoubleScalar(af);
  af = mxGetField(a, 0, "L1"); if (isDoubleRealScalar(af))  P.L1 = getDoubleScalar(af);
  af = mxGetField(a, 0, "M2"); if (isDoubleRealScalar(af))  P.m2 = getDoubleScalar(af);
  af = mxGetField(a, 0, "I2"); if (isDoubleRealScalar(af))  P.I2 = getDoubleScalar(af);
  af = mxGetField(a, 0, "L2"); if (isDoubleRealScalar(af))  P.L2 = getDoubleScalar(af);
  af = mxGetField(a, 0, "muA"); if (isDoubleRealScalar(af)) P.muA = getDoubleScalar(af);
  af = mxGetField(a, 0, "muB"); if (isDoubleRealScalar(af)) P.muB = getDoubleScalar(af);
  af = mxGetField(a, 0, "g"); if (isDoubleRealScalar(af))   P.g = getDoubleScalar(af);
  af = mxGetField(a, 0, "u"); if (isDoubleRealScalar(af))   P.u = getDoubleScalar(af);
}

void mexFunction(int nlhs, 
                 mxArray** plhs, 
                 int nrhs, 
                 const mxArray** prhs)
{
  const int arg_index_time = 0;
  const int arg_index_state = 1;
  const int arg_index_params = 2;
  const int arg_index_torque = 3;

  if ((nrhs < 3 || nrhs > 4) || nlhs != 1) {
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

  plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[arg_index_state]), mxGetN(prhs[arg_index_state]), mxREAL);
  double* ptrdz = (double *) mxGetPr(plhs[0]);

  acrobot::calculate_dotted_state(ptrdz, time, ptrz, &P);

  return;
}

double getDoubleScalar(const mxArray* a) {
  return *((const double *) mxGetPr(a));
}

bool isDoubleRealScalar(const mxArray* a) {
  return (isDoubleRealMatrix(a) && mxGetNumberOfElements(a) == 1);
}

bool isDoubleRealVector(const mxArray* a) {
  return (isDoubleRealMatrix(a) && 
          ((mxGetM(a) == 1 && mxGetN(a) >= 1) || (mxGetM(a) >= 1 && mxGetN(a) == 1)) 
          );
}

bool isDoubleRealMatrix(const mxArray* a) {
  if (a == nullptr) return false;
  if (mxIsEmpty(a)) return false;
  if (!mxIsDouble(a)) return false;
  if (mxGetNumberOfDimensions(a) != 2) return false;
  if (mxIsComplex(a)) return false;
  return true;
}
