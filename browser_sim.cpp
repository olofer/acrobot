#include <emscripten.h>
#include <cstring>
#include <cmath>
#include "acrobot.hpp"

struct AcrobotSim
{
  acrobot::params P;

  double t;

  double theta1;
  double theta2;
  double theta1dot;
  double theta2dot;

  double xy1[4];
  double xy2[4];
  double energy;

  AcrobotSim() : P(1.00, 1.50, 0.50, 1.00), 
                 t(0.0), 
                 theta1(0.0), theta2(0.0), 
                 theta1dot(0.0), theta2dot(0.0) 
  {
    compute_cms();
  }

  ~AcrobotSim()
  {}

  void compute_Is() {
    P.I1 = P.m1 * P.L1 * P.L1 / 12.0;
    P.I2 = P.m2 * P.L2 * P.L2 / 12.0;
  }

  void compute_cms() {
    const double z[4] = {theta1, theta2, theta1dot, theta2dot};
    acrobot::calculate_CM_state(xy1, xy2, z, &P);
    energy = acrobot::total_mechanical_energy(z, xy1, xy2, &P);
  }

  void evolve(double dt) {
    double k1[4];
    double k2[4];
    double k3[4];
    double k4[4];
    double z[4];

    const double t0 = t;
    const double z0[4] = {theta1, theta2, theta1dot, theta2dot};

    acrobot::calculate_dotted_state(k1, t0, z0, &P);

    const double t1 = t0 + 0.5 * dt;
    for (int i = 0; i < 4; i++)
      z[i] = z0[i] + 0.5 * dt * k1[i];

    acrobot::calculate_dotted_state(k2, t1, z, &P);

    const double t2 = t0 + 0.5 * dt;
    for (int i = 0; i < 4; i++)
      z[i] = z0[i] + 0.5 * dt * k2[i];

    acrobot::calculate_dotted_state(k3, t2, z, &P);

    const double t3 = t0 + dt;
    for (int i = 0; i < 4; i++)
      z[i] = z0[i] + dt * k3[i];

    acrobot::calculate_dotted_state(k4, t3, z, &P);

    for (int i = 0; i < 4; i++)
      z[i] = z0[i] + dt * (k1[i] / 6.0 + k2[i] / 3.0 + k3[i] / 3.0 + k4[i] / 6.0);

    t += dt;
    theta1 = z[0];
    theta2 = z[1];
    theta1dot = z[2];
    theta2dot = z[3];
    compute_cms();
  }

  double compute_freeze_torque(double y0,
                               double omega0, 
                               double zeta = 1.0) const
  {
    // Let y = (theta1 - theta2) - y0
    // Map to 2nd order linear ODE: m * y'' + c * y' + k * y = 0
    // omega0 = sqrt(k / m), c = zeta * 2 * m * omega0, zeta = 1 is critical damping
    const double y = theta1 - theta2 - y0;
    const double yprime = theta1dot - theta2dot;
    // The u-part of ypprime: ddot(theta1) - ddot(theta2) ~ (1/I1 + 1/I2) * u
    const double Ieff = (P.I1 * P.I2) / (P.I1 + P.I2);
    // y'' + (c / m) * y' + (k / m) * y = 0, can ~ be rewritten as
    // (1 / Ieff) * u + 2 * zeta * omega0 * y' + omega0^2 * y = 0; solve for u and return
    return Ieff * omega0 * (-1.0 * omega0 * y - 2.0 * zeta * yprime);
  }

  double compute_hold_torque(double y0,
                             double omega0, 
                             double zeta = 1.0) const
  {
    const double y = theta2 - y0;
    const double yprime = theta2dot;
    return P.I2 * omega0 * (2.0 * zeta * yprime + omega0 * y);
  }

};

static constexpr double nominal_gravity_acc = 9.82;

static AcrobotSim acb;


extern "C" {

EMSCRIPTEN_KEEPALIVE
void resetAcrobot()
{
  acb.t = 0.0;
  acb.theta1 = 0.0;
  acb.theta2 = 0.0;
  acb.theta1dot = 0.0;
  acb.theta2dot = 0.0;
  acb.P.m1 = 1.0;
  acb.P.L1 = 1.5;
  acb.P.m2 = 0.5;
  acb.P.L2 = 1.0;
  acb.P.g = nominal_gravity_acc;
  acb.P.u = 0.0;
  acb.P.muA = acb.P.muB = 0.0;
  acb.compute_Is();
  acb.compute_cms();
}

EMSCRIPTEN_KEEPALIVE
double getTime(void)
{
  return acb.t;
}

EMSCRIPTEN_KEEPALIVE
double getEnergy(void)
{
  return acb.energy;
}

EMSCRIPTEN_KEEPALIVE
double getTheta(int arm)
{
  return (arm == 0 ? acb.theta1 : acb.theta2);
}

EMSCRIPTEN_KEEPALIVE
double getThetaDot(int arm)
{
  return (arm == 0 ? acb.theta1dot : acb.theta2dot);
}

EMSCRIPTEN_KEEPALIVE
double getLength(int arm)
{
  return (arm == 0 ? acb.P.L1 : acb.P.L2);
}

EMSCRIPTEN_KEEPALIVE
void evolveAcrobot(double dt)
{
  acb.evolve(dt);
}

EMSCRIPTEN_KEEPALIVE
void resetAcrobotState(double th1, 
                       double th2, 
                       double th1d, 
                       double th2d)
{
  acb.t = 0.0;
  acb.theta1 = th1;
  acb.theta2 = th2;
  acb.theta1dot = th1d;
  acb.theta2dot = th2d;
  acb.compute_cms();
  acb.P.u = 0.0;
}

EMSCRIPTEN_KEEPALIVE
void applyTorque(double tau)
{
  acb.P.u = tau;
}

EMSCRIPTEN_KEEPALIVE
double getTorque(void)
{
  return acb.P.u;
}

EMSCRIPTEN_KEEPALIVE
double freezeTorque(double y0, 
                    double omega0, 
                    double zeta)
{
  return acb.compute_freeze_torque(y0, omega0, zeta);
}

EMSCRIPTEN_KEEPALIVE
double holdTorque(double y0, 
                  double omega0, 
                  double zeta)
{
  return acb.compute_hold_torque(y0, omega0, zeta);
}

EMSCRIPTEN_KEEPALIVE
void increaseFriction(int joint, double inc)
{
  if (joint == 0) {
    acb.P.muA += inc;
    if (acb.P.muA < 0.0) acb.P.muA = 0.0;
  }
  else {
    acb.P.muB += inc;
    if (acb.P.muB < 0.0) acb.P.muB = 0.0;
  }
}

EMSCRIPTEN_KEEPALIVE
double getFriction(int joint)
{
  return (joint == 0 ? acb.P.muA : acb.P.muB);
}

EMSCRIPTEN_KEEPALIVE
void setFriction(int joint, double mu)
{
  if (joint == 0) {
    acb.P.muA = mu;
  } else if (joint == 1) {
    acb.P.muB = mu;
  }
}

EMSCRIPTEN_KEEPALIVE
double cosineAlpha(double theta1, double theta2)
{
  double CoM[2];
  double cosa = 0.0;
  acrobot::total_cartesian_CoM(CoM, &cosa, theta1, theta2, &acb.P);
  return cosa;
}

EMSCRIPTEN_KEEPALIVE
double directionIndicator(double alpha) // +1 if climbing up, -1 if falling down (or still)
{
  const double ydot = std::cos(acb.theta1 + alpha) * acb.theta1dot;
  return (ydot > 0.0 ? +1.0 : -1.0);
}

} // close extern "C"
