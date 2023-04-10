/*
 * USAGE: USAGE: [V, A, g1, g2, g3, g4] = acrobot_dpsolve(P, npts, itrs, dt, tdisc);
 *
 */ 

#include "mex.h"
#include <cstring>
#include <cmath>

#include "../acrobot.hpp"
#include "acrobot_mex_utils.hpp"

#include <vector>
#include <limits>

#include <omp.h>

const double _one_pi = 3.14159265358979323846;

template <typename T>
class acrobotDP
{
private:
  std::vector<double> grid_th1;
  std::vector<double> grid_th2;
  std::vector<double> grid_th1d;
  std::vector<double> grid_th2d;

  int dims[4];
  double deltas[4];

  std::vector<T> value;
  std::vector<T> value_update;
  std::vector<int8_t> action;
  double time;

  T edge_value;

public:
  acrobotDP(int nth1, 
            int nth2, 
            int ndot1, 
            int ndot2,
            double maxdot1 = 4.0 * _one_pi,
            double maxdot2 = 4.0 * _one_pi) : 
  grid_th1(nth1), 
  grid_th2(nth2), 
  grid_th1d(ndot1), 
  grid_th2d(ndot2),
  value(nth1 * nth2 * ndot1 * ndot2), 
  value_update(value),
  action(nth1 * nth2 * ndot1 * ndot2), 
  dims{nth1, nth2, ndot1, ndot2},
  time(0.0),
  edge_value(-100.0)
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

  void set_edge_value(T ev) {
    edge_value = ev;
  }

  int lookupindex(int i0, int i1, int i2, int i3) const {
    if (i2 < 0 || i2 >= dims[2]) return edge_value;
    if (i3 < 0 || i3 >= dims[3]) return edge_value;

    if (i0 < 0) i0 += dims[0];
    if (i0 >= dims[0]) i0 -= dims[0];

    if (i1 < 0) i1 += dims[1];
    if (i1 >= dims[1]) i1 -= dims[1];

    /*if (i2 < 0) i2 = 0;
    if (i2 >= dims[2]) i2 = dims[2] - 1;

    if (i3 < 0) i3 = 0;
    if (i3 >= dims[3]) i3 = dims[3] - 1; */

    return sub2ind(i0, i1, i2, i3);
  }

  T lookup(int i0, int i1, int i2, int i3) const {
    return value[lookupindex(i0, i1, i2, i3)];
  }

  void assign(int i0, int i1, int i2, int i3, T q) {
    value[lookupindex(i0, i1, i2, i3)] = q;
  }

  void impose(int i0, int i1, int i2, int i3, T q) {
    value[lookupindex(i0, i1, i2, i3)] += q;
  }

  void griddify(const double* x,
                int *i,
                double* eta) const
  {
    const double y0 = (x[0] - grid_th1[0]) / deltas[0];
    i[0] = (int) std::floor(y0);
    eta[0] = y0 - i[0];

    const double y1 = (x[1] - grid_th2[0]) / deltas[1];
    i[1] = (int) std::floor(y1);
    eta[1] = y1 - i[1];

    const double y2 = (x[2] - grid_th1d[0]) / deltas[2];
    i[2] = (int) std::floor(y2);
    eta[2] = y2 - i[2];

    const double y3 = (x[3] - grid_th2d[0]) / deltas[3];
    i[3] = (int) std::floor(y3);
    eta[3] = y3 - i[3];
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
    int i[4];
    double eta[4];
    griddify(x, i, eta);
    return interp4d(i[0], eta[0], 
                    i[1], eta[1],
                    i[2], eta[2],
                    i[3], eta[3]);
  }

  void scatter(int i0, double eta0, // assume 0 <= eta < 1
               int i1, double eta1,
               int i2, double eta2,
               int i3, double eta3,
               T val)
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
               T mass)
  { 
    int i[4];
    double eta[4];
    griddify(x, i, eta);
    scatter(i[0], eta[0], 
            i[1], eta[1],
            i[2], eta[2],
            i[3], eta[3],
            mass);
  }

  double update_time() const { return time; }

  bool update(const acrobot::params* P,
              const std::vector<double>& ulevels,
              double dt,
              double tdisc,
              bool use_euler = false) 
  {
    acrobot::params localP(*P);

    const double gamma = std::exp(-1.0 * dt / tdisc); // 1-step (dt) survival probability; rate = 1/tdisc

    const double epsa = 1.0e-12;
    const double lowest_value = std::numeric_limits<double>::lowest();
    const int numlevels = ulevels.size();

    double updatesum = 0.0;
    double sum_value_delta = 0.0;
    int action_edits = 0;

    #pragma omp parallel for firstprivate(localP) reduction(+:updatesum,sum_value_delta,action_edits)
    for (int k = 0; k < size(); k++) {
      int ik[4];
      ind2sub(k, ik);
      const double xk[4] = {grid_th1[ik[0]], 
                            grid_th2[ik[1]], 
                            grid_th1d[ik[2]], 
                            grid_th2d[ik[3]]};

      const double Rk = reward_swing(xk, &localP, 0.0);  // assume independence of u(a)

      double max_value = lowest_value;
      double max_value_clean = lowest_value;
      int argmax = -1;

      for (int a = 0; a < numlevels; a++) {
        double xnext[4];
        const double ua = ulevels[a];
        localP.u = ua;

        if (use_euler) {
          acrobot::step_euler(xnext, xk, dt, &localP);
        } else {
          acrobot::step_heun(xnext, xk, dt, &localP);
        }

        const double vnext_a_node = dt * Rk + gamma * interp4d(xnext); 
        const double vnext_a_regu = vnext_a_node - epsa * ua * ua;

        if (vnext_a_regu > max_value) {
          max_value = vnext_a_regu;
          max_value_clean = vnext_a_node;
          argmax = a;
        }
      }

      value_update[k] = static_cast<T>(max_value_clean);
      sum_value_delta += std::fabs(value_update[k] - value[k]);

      if (action[k] != (int8_t) argmax) {
        action[k] = (int8_t) argmax;
        action_edits++;
      }

      updatesum += max_value_clean;
    }

    overwrite();

    time += dt;

    mexPrintf("[%s]: <value>=%e (lvls=%i, aedits=%i); bkwdtm=%f, eps(delta)=%e\n", 
              __func__, updatesum / size(), numlevels, action_edits, time, sum_value_delta / std::fabs(updatesum));
    
    return true;
  }

  double wrap_angle_diff(double a) const {
    while (a < -_one_pi) a += 2.0 * _one_pi;
    while (a > _one_pi) a -= 2.0 * _one_pi;
    return a;
  }

/*
  void initialize_terminal(int cfgnum) 
  {
    const double equilibria[4 * 4] = {-_one_pi / 2.0, -_one_pi / 2.0, 0.0, 0.0,
                                      -_one_pi / 2.0,  _one_pi / 2.0, 0.0, 0.0,
                                       _one_pi / 2.0, -_one_pi / 2.0, 0.0, 0.0,
                                       _one_pi / 2.0,  _one_pi / 2.0, 0.0, 0.0};
    double target[4];
    if (cfgnum < 0 || cfgnum > 3) cfgnum = 0;
    std::memcpy(target, &equilibria[4 * cfgnum], 4 * sizeof(double));
    for (int k = 0; k < size(); k++) {
      int ik[4];
      ind2sub(k, ik);
      const double xk[4] = {grid_th1[ik[0]], grid_th2[ik[1]], grid_th1d[ik[2]], grid_th2d[ik[3]]};

      const double delta0 = wrap_angle_diff(target[0] - xk[0]);
      const double delta1 = wrap_angle_diff(target[1] - xk[1]);
      const double delta2 = target[2] - xk[2];
      const double delta3 = target[3] - xk[3];

      double reward = edge_value;

      if (std::fabs(delta0) < 2.0 * _one_pi / 12 &&
          std::fabs(delta1) < 2.0 * _one_pi / 12 &&
          std::fabs(delta2) < 0.50 &&
          std::fabs(delta3) < 0.50)
        reward = 0.0;

      value[k] = reward;
    }
  }
*/

  double reward_swing(const double* z, 
                      const acrobot::params* P, 
                      double E) const
  {
    double xy1[4];
    double xy2[4];

    acrobot::calculate_CM_state(xy1, xy2, z, P);
    const double Ek = acrobot::total_mechanical_energy(z, xy1, xy2, P);
    const double Esize = 1.0 + std::fabs(E);

    const double v1x = std::cos(z[0]);
    const double v1y = std::sin(z[0]);
    const double v2x = std::cos(z[1]);
    const double v2y = std::sin(z[1]);

    const double dot12 = v1x * v2x + v1y * v2y; // |.| <= 1, want it to be 1

    return -1.0 * (std::fabs(Ek - E) / Esize + (1.0 - dot12));
  }

  void clear() {
    time = 0.0;
    std::memset(value.data(), 0, sizeof(T) * value.size());
    std::memset(value_update.data(), 0, sizeof(T) * value_update.size());
    std::memset(action.data(), 0, sizeof(int8_t) * action.size());
  }

  void test_scatter_interp(double putval) {

    const double X0[4] = {_one_pi / 2.0, _one_pi / 2.0, 0.0, 0.0};
    const double X1[4] = {-_one_pi / 2.0, -_one_pi / 2.0, 0.0, 0.0};
    const double X2[4] = {0.0, 0.0, 0.0, 0.0};
    const double X3[4] = {0.0, 0.0, 1.0, -1.0};
    const double X4[4] = {1.11, 2.22, 3.33, 4.44};

    const double* X[5] = {X0, X1, X2, X3, X4};

    const double eta[4] = {0.5, 0.12345, 0.2523, 0.96};

    for (int j = 0; j < 5; j++) {
      clear();
      scatter(X[j], putval);
      const double getval = interp4d(X[j]);
      const double sj = totalValueSum();
      const double sjo = totalValueSumOffsetGrid(eta);
      mexPrintf("[%s]: put value = %f (%i)\n", __func__, putval, j);
      mexPrintf("[%s]: total array sum = %f, offset sum = %f (%i)\n", __func__, sj, sjo, j);
      mexPrintf("[%s]: get value = %f (%i)\n", __func__, getval, j);
      mexPrintf("--- --- ---\n");
    }
  }

  int sizeofValueType() const { return sizeof(T); }

  T valueAt(int k) const { return value[k]; }

  int8_t actionAt(int k) const { return action[k]; }

  double totalValueSum() const {
    double s = 0.0;
    for (int k = 0; k < size(); k++) {
      s += value[k];
    }
    return s;
  }

  double totalValueSumOffsetGrid(const double* eta) const {
    int ik[4];
    double s = 0.0;
    for (int k = 0; k < size(); k++) {
      ind2sub(k, ik);
      const double Vk = interp4d(ik[0], eta[0], 
                                 ik[1], eta[1],
                                 ik[2], eta[2],
                                 ik[3], eta[3]);
      s += Vk;
    }
    return s;
  }

  double gridPoint(int d, int i) const {
    switch (d) {
      case 0:
        return grid_th1[i];
      case 1:
        return grid_th2[i];
      case 2:
        return grid_th1d[i];
      case 3:
        return grid_th2d[i];
      default:
        return std::numeric_limits<double>::quiet_NaN();
    }
  }

  double averageCellDistance(const acrobot::params* P,
                             double dt) const
  {
    acrobot::params localP(*P);
    localP.u = 0.0;
    double total_deltas[4] = {0.0, 0.0, 0.0, 0.0};
    for (int k = 0; k < size(); k++) {
      int ik[4];
      ind2sub(k, ik);
      const double xk[4] = { grid_th1[ik[0]], grid_th2[ik[1]], grid_th1d[ik[2]], grid_th2d[ik[3]] };
      double xnext[4];
      acrobot::step_heun(xnext, xk, dt, &localP);
      for (int d = 0; d < 4; d++) {
        total_deltas[d] += std::fabs(xnext[d] - xk[d]);
      }
    }
    double avg_abs_deltas[4] = { total_deltas[0] / size(),
                                 total_deltas[1] / size(),
                                 total_deltas[2] / size(),
                                 total_deltas[3] / size() };
    return ( avg_abs_deltas[0] / deltas[0] + 
             avg_abs_deltas[1] / deltas[1] +
             avg_abs_deltas[2] / deltas[2] +
             avg_abs_deltas[3] / deltas[3] ) / 4;
  }

private:

  void overwrite() {
    std::memcpy(value.data(), value_update.data(), sizeof(T) * value.size());
  }

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
  const int arg_index_disctime = 4;

  if (nrhs != 5 || nlhs > 6) {
    mexErrMsgTxt("USAGE: [V, A, g1, g2, g3, g4] = acrobot_dpsolve(P, npts, itrs, dt, tdisc);");
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

  if (!isDoubleRealScalar(prhs[arg_index_disctime])) {
    mexErrMsgTxt("Expects tdisc to be real scalar");
  }

  const int itrs = (int) std::round(getDoubleScalar(prhs[arg_index_iterations]));

  if (itrs < 0) {
    mexErrMsgTxt("itrs >= 0 required");
  }

  const double dt = getDoubleScalar(prhs[arg_index_deltatime]);

  if (dt <= 0.0) {
    mexErrMsgTxt("dt > 0 required");
  }

  const double tdisc = getDoubleScalar(prhs[arg_index_disctime]);
  
  if (tdisc <= 0.0) {
    mexErrMsgTxt("tdisc > 0 required");
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

  //acrobotDP<double> acbdp(npts[0], npts[1], npts[2], npts[3]);
  acrobotDP<float> acbdp(npts[0], npts[1], npts[2], npts[3]);

  mexPrintf("num. cells = %i\n", acbdp.size());
  mexPrintf("dims       = %i, %i, %i, %i\n", acbdp.dim(0), acbdp.dim(1), acbdp.dim(2), acbdp.dim(3));
  mexPrintf("deltas     = %f, %f, %f, %f\n", acbdp.delta(0), acbdp.delta(1), acbdp.delta(2), acbdp.delta(3));
  mexPrintf("sizeof(T)  = %i (sizeof value element type)\n", acbdp.sizeofValueType());
  mexPrintf("dt, tdisc  = %f [sec], %f [sec]\n", dt, tdisc);
  mexPrintf("dt * itrs  = %f [sec] (itrs = %i)\n", dt * itrs, itrs);

  if (!acbdp.check_ind2sub2ind()) {
    mexErrMsgTxt("Self-check of ind2sub, sub2ind failed");
  }

  if (nlhs != 6) {
    mexPrintf("[%s]: (not 6 output arguments) running self-check\n", __func__);
    acbdp.test_scatter_interp(2.56);

    const double avg_cell_dist = acbdp.averageCellDistance(&P, dt);
    mexPrintf("[%s]: <cell-delta> = %f\n", __func__, avg_cell_dist);

    return;
  }

  acbdp.clear();
  acbdp.set_edge_value(-2.0 * itrs * dt);

  std::vector<double> ulevels = {-1.0, 0.0, +1.0};

  for (int i = 0; i < itrs; i++) {
    if (!acbdp.update(&P, ulevels, dt, tdisc)) break;
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

  for (int i = 0; i < 4; i++) {
    plhs[2 + i] = mxCreateDoubleMatrix(acbdp.dim(i), 1, mxREAL);
    double* outi = (double *) mxGetPr(plhs[2 + i]);
    for (int j = 0; j < acbdp.dim(i); j++) {
      outi[j] = acbdp.gridPoint(i, j);
    }
  }

  return;
}
