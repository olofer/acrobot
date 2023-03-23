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
#include <limits>

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
    deltas[0] = 2.0 * _one_pi / nth1;
    deltas[1] = 2.0 * _one_pi / nth2;
    deltas[2] = 2.0 * maxdot1 / (ndot1 - 1);
    deltas[3] = 2.0 * maxdot2 / (ndot2 - 1);
  }

  ~acrobotDP() {}

  int dim(int i) const { return dims[i]; }

  double delta(int i) const { return deltas[i]; }

  int size() const { return value.size(); }

  int sub2ind(int i0, int i1, int i2, int i3) const {
    return i0 + dims[0] * (i1 + dims[1] * (i2 + dims[2] * i3));
  }

  int sub2ind(const int* i) const {
    return sub2ind(i[0], i[1], i[2], i[3]); 
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

  double lookup(int i0, int i1, int i2, int i3) const {
    if (i2 < 0 || i2 >= dims[2]) return 0.0;
    if (i3 < 0 || i3 >= dims[3]) return 0.0;
    if (i0 < 0) i0 += dims[0];
    if (i0 >= dims[0]) i0 -= dims[0];
    if (i1 < 0) i1 += dims[1];
    if (i1 >= dims[1]) i1 -= dims[1];
    return value[sub2ind(i0, i1, i2, i3)];
  }

  void interp4d_weights(const double* eta, double* w) const {
    const double w0 = (1.0 - eta[0]);
    const double w1 = (eta[0]);

    const double w00 = w0 * (1.0 - eta[1]);
    const double w01 = w0 * (eta[1]);
    const double w10 = w1 * (1.0 - eta[1]);
    const double w11 = w1 * (eta[1]);

    const double w000 = w00 * (1.0 - eta[2]);
    const double w001 = w00 * (eta[2]);
    const double w010 = w01 * (1.0 - eta[2]);
    const double w011 = w01 * (eta[2]);
    const double w100 = w10 * (1.0 - eta[2]);
    const double w101 = w10 * (eta[2]);
    const double w110 = w11 * (1.0 - eta[2]);
    const double w111 = w11 * (eta[2]);

    w[0]  = w000 * (1.0 - eta[3]);  // w{0,0,0,0}
    w[1]  = w000 * (eta[3]);        // w{0,0,0,1}
    w[2]  = w001 * (1.0 - eta[3]);  // ...
    w[3]  = w001 * (eta[3]);
    w[4]  = w010 * (1.0 - eta[3]);
    w[5]  = w010 * (eta[3]);
    w[6]  = w011 * (1.0 - eta[3]);
    w[7]  = w011 * (eta[3]);

    w[8]  = w100 * (1.0 - eta[3]);
    w[9]  = w100 * (eta[3]);
    w[10] = w101 * (1.0 - eta[3]);
    w[11] = w101 * (eta[3]);
    w[12] = w110 * (1.0 - eta[3]);
    w[13] = w110 * (eta[3]);        // ...
    w[14] = w111 * (1.0 - eta[3]);  // w{1,1,1,0}
    w[15] = w111 * (eta[3]);        // w{1,1,1,1}
  }

  double interp4d(int i0, double eta0, // assume 0 <= eta < 1
                  int i1, double eta1,
                  int i2, double eta2,
                  int i3, double eta3) const {
    
    const double eta[4] = {eta0, eta1, eta2, eta3};

    double weights[16];
    interp4d_weights(eta, weights);

    double values[16] = {
      lookup(i0, i1, i2, i3),
      lookup(i0, i1, i2, i3 + 1),
      lookup(i0, i1, i2 + 1, i3),
      lookup(i0, i1, i2 + 1, i3 + 1),

      lookup(i0, i1 + 1, i2, i3),
      lookup(i0, i1 + 1, i2, i3 + 1),
      lookup(i0, i1 + 1, i2 + 1, i3),
      lookup(i0, i1 + 1, i2 + 1, i3 + 1),

      lookup(i0 + 1, i1, i2, i3),
      lookup(i0 + 1, i1, i2, i3 + 1),
      lookup(i0 + 1, i1, i2 + 1, i3),
      lookup(i0 + 1, i1, i2 + 1, i3 + 1),

      lookup(i0 + 1, i1 + 1, i2, i3),
      lookup(i0 + 1, i1 + 1, i2, i3 + 1),
      lookup(i0 + 1, i1 + 1, i2 + 1, i3),
      lookup(i0 + 1, i1 + 1, i2 + 1, i3 + 1)
    };

    double sum = 0.0;
    for (int i = 0; i < 16; i++) {
      sum += weights[i] * values[i];
    }

    return sum;
  }

  double interp4d(const double* x) const {
    const double y0 = (x[0] - grid_th1[0]) / deltas[0];
    const double i0 = std::floor(y0);

    const double y1 = (x[1] - grid_th2[0]) / deltas[1];
    const double i1 = std::floor(y1);

    const double y2 = (x[2] - grid_th1d[0]) / deltas[2];
    const double i2 = std::floor(y2);

    const double y3 = (x[3] - grid_th2d[0]) / deltas[3];
    const double i3 = std::floor(y3);

    return interp4d((int) i0, y0 - i0, 
                    (int) i1, y1 - i1,
                    (int) i2, y2 - i2,
                    (int) i3, y3 - i3);
  }

  // FIXME: "add(x, weight)" function; add weight to value cells centered at location x (~ reverse of interp4d)

  void euler_update(const acrobot::params* P,
                    const std::vector<double>& ulevels,
                    double dt) const 
  {
    acrobot::params localP(*P);
    int ik[4];
    double xdot[4];

    if (localP.L2 != P->L2 || localP.I1 != P->I1) return;

    const int numlevels = ulevels.size();
    double totalsum = 0.0;

    for (int k = 0; k < size(); k++) {
      ind2sub(k, ik);
      const double xk[4] = {grid_th1[ik[0]], 
                            grid_th2[ik[1]], 
                            grid_th1d[ik[2]], 
                            grid_th2d[ik[3]]};

      double max_value = std::numeric_limits<double>::lowest();
      int argmax = -1;

      for (int a = 0; a < numlevels; a++) {
        localP.u = ulevels[a];
        acrobot::calculate_dotted_state(xdot, 0.0, xk, &localP);

        const double xnext[4] = {xk[0] + dt * xdot[0], 
                                 xk[1] + dt * xdot[1], 
                                 xk[2] + dt * xdot[2], 
                                 xk[3] + dt * xdot[3]};

        const double vnext_a = interp4d(xnext);

        if (vnext_a > max_value) {
          max_value = vnext_a;
          argmax = a;
        } // TODO: on equality; select a if ulevels[a] < current selection (if any, i.e. argmax != -1)

        totalsum += vnext_a;
      }

      // Here we need to select the best value and assign to action[.]
      // store argmax or ulevels[argmax] as action
    }

    mexPrintf("[%s]: sum(all next values) = %e (%i actions)\n", __func__, totalsum, numlevels);
  }

  int initialize() {
    std::memset(value.data(), 0, sizeof(float) * value.size());
    // FIXME: assign value one to target cell close to: {pi/2,pi/2,0,0}
    value[sub2ind(10, 11, 12, 13)] = 1.0;
    return 0;
  }

  int update(const acrobot::params* P, 
             double dt) 
  {
    double sum = 0.0;
    int ik[4];

    for (int k = 0; k < size(); k++) {
      ind2sub(k, ik);
      const double Vk = interp4d(ik[0], 0.125, 
                                 ik[1], 0.25,
                                 ik[2], 0.50,
                                 ik[3], 0.75);
      sum += Vk;
      if (Vk != 0.0) {
        mexPrintf("Vk = %e @ k = %i\n", Vk, k);
      }
    }
    mexPrintf("[%s]: sum = %e\n", __func__, sum);
    return 0;
  }

  // FIXME: missing member functions to "apply" the control law: i.e. control = evaluate(state).

private:
  std::vector<double> grid_th1;
  std::vector<double> grid_th2;
  std::vector<double> grid_th1d;
  std::vector<double> grid_th2d;

  int dims[4];
  double deltas[4];

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
  mexPrintf("dims       = %i, %i, %i, %i\n", acbdp.dim(0), acbdp.dim(1), acbdp.dim(2), acbdp.dim(3));
  mexPrintf("deltas     = %f, %f, %f, %f\n", acbdp.delta(0), acbdp.delta(1), acbdp.delta(2), acbdp.delta(3));

  if (!acbdp.check_ind2sub2ind()) {
    mexErrMsgTxt("Self-check of ind2sub, sub2ind failed");
  }
  
  acbdp.initialize();
  acbdp.update(&P, 1.0e-3);

  std::vector<double> ulevels = {-1.0, 0.0, +1.0};
  acbdp.euler_update(&P, ulevels, 1.0e-3);

  const mwSize dims[4] = {acbdp.dim(0), acbdp.dim(1), acbdp.dim(2), acbdp.dim(3)};
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
