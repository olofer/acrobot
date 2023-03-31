#pragma once

#if defined(HAVE_OCTAVE) || defined(MATLAB_MEX_FILE)

bool isDoubleRealMatrix(const mxArray* a) {
  if (a == nullptr) return false;
  if (mxIsEmpty(a)) return false;
  if (!mxIsDouble(a)) return false;
  if (mxGetNumberOfDimensions(a) != 2) return false;
  if (mxIsComplex(a)) return false;
  return true;
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

#endif
