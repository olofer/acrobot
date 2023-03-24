/*
 * USAGE: [V, A] = acrobot_dpsolve(P, npts, itrs, dt);
 *
 */

// FIXME: return grid vectors
// (and update test script; including visualization; and apply feedback law)
// FIXME: implement Heuns integration method as alternative; 
// FIXME: must be able to specify the min/max solver bounds for the ang. velocities
// FIXME: determine the operational max(abs(thetadot)) from the JS/WASM simulation
// ..

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
private:
  std::vector<double> grid_th1;
  std::vector<double> grid_th2;
  std::vector<double> grid_th1d;
  std::vector<double> grid_th2d;

  int dims[4];
  double deltas[4];

  std::vector<float> value;
  std::vector<float> value_update;
  std::vector<int8_t> action;
  double time;

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
  value_update(value),
  action(nth1 * nth2 * ndot1 * ndot2), 
  dims{nth1, nth2, ndot1, ndot2},
  time(0.0)
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

  int lookupindex(int i0, int i1, int i2, int i3) const {
    if (i2 < 0 || i2 >= dims[2]) return 0.0;
    if (i3 < 0 || i3 >= dims[3]) return 0.0;
    if (i0 < 0) i0 += dims[0];
    if (i0 >= dims[0]) i0 -= dims[0];
    if (i1 < 0) i1 += dims[1];
    if (i1 >= dims[1]) i1 -= dims[1];
    return sub2ind(i0, i1, i2, i3);
  }

  double lookup(int i0, int i1, int i2, int i3) const {
    return value[lookupindex(i0, i1, i2, i3)];
  }

  void assign(int i0, int i1, int i2, int i3, double q) {
    value[lookupindex(i0, i1, i2, i3)] = q;
  }

  void impose(int i0, int i1, int i2, int i3, double q) {
    value[lookupindex(i0, i1, i2, i3)] += q;
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

  void scatter(int i0, double eta0, // assume 0 <= eta < 1
               int i1, double eta1,
               int i2, double eta2,
               int i3, double eta3,
               double val)
  {
    const double eta[4] = {eta0, eta1, eta2, eta3};

    double weights[16];
    interp4d_weights(eta, weights);

    impose(i0, i1, i2, i3, val * weights[0]);
    impose(i0, i1, i2, i3 + 1, val * weights[1]);
    impose(i0, i1, i2 + 1, i3, val * weights[2]);
    impose(i0, i1, i2 + 1, i3 + 1, val * weights[3]);

    impose(i0, i1 + 1, i2, i3, val * weights[4]);
    impose(i0, i1 + 1, i2, i3 + 1, val * weights[5]);
    impose(i0, i1 + 1, i2 + 1, i3, val * weights[6]);
    impose(i0, i1 + 1, i2 + 1, i3 + 1, val * weights[7]);

    impose(i0 + 1, i1, i2, i3, val * weights[8]);
    impose(i0 + 1, i1, i2, i3 + 1, val * weights[9]);
    impose(i0 + 1, i1, i2 + 1, i3, val * weights[10]);
    impose(i0 + 1, i1, i2 + 1, i3 + 1, val * weights[11]);

    impose(i0 + 1, i1 + 1, i2, i3, val * weights[12]);
    impose(i0 + 1, i1 + 1, i2, i3 + 1, val * weights[13]);
    impose(i0 + 1, i1 + 1, i2 + 1, i3, val * weights[14]);
    impose(i0 + 1, i1 + 1, i2 + 1, i3 + 1, val * weights[15]);
  }

  void scatter(const double* x, 
               double mass)
  {
    const double y0 = (x[0] - grid_th1[0]) / deltas[0];
    const double i0 = std::floor(y0);

    const double y1 = (x[1] - grid_th2[0]) / deltas[1];
    const double i1 = std::floor(y1);

    const double y2 = (x[2] - grid_th1d[0]) / deltas[2];
    const double i2 = std::floor(y2);

    const double y3 = (x[3] - grid_th2d[0]) / deltas[3];
    const double i3 = std::floor(y3);

    scatter((int) i0, y0 - i0, 
            (int) i1, y1 - i1,
            (int) i2, y2 - i2,
            (int) i3, y3 - i3,
            mass);
  }

  double update_time() const { return time; }

  void euler_update(const acrobot::params* P,
                    const std::vector<double>& ulevels,
                    double dt) 
  {
    acrobot::params localP(*P);
    int ik[4];
    double xdot[4];
    int actionChanges = 0;

    const int numlevels = ulevels.size();
    double totalsum = 0.0;
    double updatesum = 0.0;

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
        } else if (vnext_a == max_value && ulevels[a] == 0.0) {
          argmax = a;
        }

        totalsum += vnext_a;
      }

      value_update[k] = max_value;
      if (action[k] != (int8_t) argmax) {
        action[k] = (int8_t) argmax;
        actionChanges++;
      }

      updatesum += max_value;

      // Here we need to select the best value and assign to action[.]
      // store argmax or ulevels[argmax] as action
    }

    std::memcpy(value.data(), value_update.data(), sizeof(float) * value.size());
    time += dt;

    mexPrintf("[%s]: totalsum=%e (across %i actions); sum(update)=%e (changes=%i); bkwdtm=%f\n", __func__, totalsum, numlevels, updatesum, actionChanges, time);
  }

  int initialize() {
    time = 0.0;
    std::memset(value.data(), 0, sizeof(float) * value.size());
    std::memset(value_update.data(), 0, sizeof(float) * value_update.size());
    std::memset(action.data(), 0, sizeof(int8_t) * action.size());

    const double xupright[4] = {_one_pi / 2.0, _one_pi / 2.0, 0.0, 0.0};
    scatter(xupright, 1.0);

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

  float valueAt(int k) const { return value[k]; }

  int8_t actionAt(int k) const { return action[k]; }
};

void mexFunction(int nlhs, 
                 mxArray** plhs, 
                 int nrhs, 
                 const mxArray** prhs)
{
  const int arg_index_params = 0;
  const int arg_index_gridpoints = 1;
  const int arg_index_iterations = 2;
  const int arg_index_deltatime = 3;

  if (nrhs != 4 || nlhs != 2) {
    mexErrMsgTxt("USAGE: [V, A] = acrobot_dpsolve(P, npts, itrs, dt);");
  }

  if (!(isDoubleRealVector(prhs[arg_index_gridpoints]) && 
        mxGetNumberOfElements(prhs[arg_index_gridpoints]) == 4)) {
    mexErrMsgTxt("npts should be a 4d vector");
  }

  if (!isDoubleRealScalar(prhs[arg_index_iterations])) {
    mexErrMsgTxt("Expects itrs to be real scalar");
  }

  if (!isDoubleRealScalar(prhs[arg_index_deltatime])) {
    mexErrMsgTxt("Expects dt to be real scalar");
  }

  const int itrs = (int) std::round(getDoubleScalar(prhs[arg_index_iterations]));

  if (itrs < 0) {
    mexErrMsgTxt("itrs >= 0 required");
  }

  const double dt = getDoubleScalar(prhs[arg_index_deltatime]);
  if (dt <= 0.0) {
    mexErrMsgTxt("dt > 0 required");
  }

  const double* pdims = mxGetPr(prhs[arg_index_gridpoints]);
  const int npts[4] = {(int) std::floor(pdims[0]),
                       (int) std::floor(pdims[1]),
                       (int) std::floor(pdims[2]),
                       (int) std::floor(pdims[3])};

  for (int i = 0; i < 4; i++) {
    if (npts[i] < 3) {
      mexErrMsgTxt("At least 3 points per dimension is expected");
    }
  }

  if (!mxIsStruct(prhs[arg_index_params])) {
    mexErrMsgTxt("Expects P to be a (scalar) struct");
  }

  acrobot::params P;
  set_from_struct(P, prhs[arg_index_params]);
  if (!P.sanityChecksOut()) {
    mexErrMsgTxt("P does not hold a complete and valid set of parameters");
  }

  //acrobotDP acbdp(36, 34, 33, 35);
  acrobotDP acbdp(npts[0], npts[1], npts[2], npts[3]);

  mexPrintf("num. cells = %i\n", acbdp.size());
  mexPrintf("dims       = %i, %i, %i, %i\n", acbdp.dim(0), acbdp.dim(1), acbdp.dim(2), acbdp.dim(3));
  mexPrintf("deltas     = %f, %f, %f, %f\n", acbdp.delta(0), acbdp.delta(1), acbdp.delta(2), acbdp.delta(3));

  if (!acbdp.check_ind2sub2ind()) {
    mexErrMsgTxt("Self-check of ind2sub, sub2ind failed");
  }
  
  acbdp.initialize();
  acbdp.update(&P, dt);

  std::vector<double> ulevels = {-1.0, -0.5, 0.0, 0.5, +1.0};
  for (int i = 0; i < itrs; i++) {
    acbdp.euler_update(&P, ulevels, dt);
  }

  // Output 4D arrays have same memory layout as Octave so just linear copy
  const mwSize dims[4] = {acbdp.dim(0), acbdp.dim(1), acbdp.dim(2), acbdp.dim(3)};
  plhs[0] = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);

  double* V = (double*) mxGetPr(plhs[0]);
  double* A = (double*) mxGetPr(plhs[1]);

  for (int k = 0; k < acbdp.size(); k++) {
    V[k] = acbdp.valueAt(k);
    A[k] = ulevels[acbdp.actionAt(k)];
  }

  /*for (int k = 0; k < acbdp.size(); k++) {
    int indices[4];
    acbdp.ind2sub(k, indices);
    V[acbdp.sub2ind(indices[0], indices[1], indices[2], indices[3])] = (double) k;
    A[acbdp.sub2ind(indices[0], indices[1], indices[2], indices[3])] = (double) (acbdp.size() - k - 1);
  }*/

  return;
}
