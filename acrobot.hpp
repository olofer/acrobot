#pragma once

namespace acrobot
{

struct params {
  double m1;
  double I1;
  double L1;
  double m2;
  double I2;
  double L2;
  double muA;
  double muB;
  double g;
  double u;

  params() {
    std::memset(this, 0, sizeof(params));
  }

  params(double m1, 
         double L1, 
         double m2, 
         double L2) : m1(m1), L1(L1), m2(m2), L2(L2)
  {
    I1 = m1 * L1 * L1 / 12.0; // thin rod moment of inertia
    I2 = m2 * L2 * L2 / 12.0;
    muA = muB = 0.0;
    g = 9.82;
    u = 0.0;
  }

  bool sanityChecksOut() const {
    return (m1 > 0 && I1 > 0 && L1 > 0 && 
            m2 > 0 && I2 > 0 && L2 > 0 && 
            muA >= 0 && muB >= 0);
  }
};

void pivotfree_solve_4x4(const double* a, 
                         const double* b,
                         double* x)
{
  // Factorize A = L*U with upper-triangular U (with ones on diagonal) and lower-triangular L
  const double l11 = a[0];
  const double u12 = a[1] / l11;
  const double u13 = a[2] / l11;
  const double u14 = a[3] / l11;
  const double l21 = a[4];
  const double l22 = a[5] - l21 * u12;
  const double u23 = (a[6] - l21 * u13) / l22;
  const double u24 = (a[7] - l21 * u14) / l22;
  const double l31 = a[8];
  const double l32 = a[9] - l31 * u12;
  const double l33 = a[10] - l31 * u13 - l32 * u23;
  const double u34 = (a[11] - l31 * u14 - l32 * u24) / l33;
  const double l41 = a[12];
  const double l42 = a[13] - l41 * u12;
  const double l43 = a[14] - l41 * u13 - l42 * u23;
  const double l44 = a[15] - l41 * u14 - l42 * u24 - l43 * u34;
  // Backsubstitutions: b = L*y, y = U*x
  const double y1 = (b[0]) / l11;
  const double y2 = (b[1] - l21 * y1) / l22;
  const double y3 = (b[2] - l31 * y1 - l32 * y2) / l33;
  const double y4 = (b[3] - l41 * y1 - l42 * y2 - l43 * y3) / l44;
  x[3] = y4;
  x[2] = y3 - u34 * x[3];
  x[1] = y2 - u23 * x[2] - u24 * x[3];
  x[0] = y1 - u12 * x[1] - u13 * x[2] - u14 * x[3];
}

void calculate_dotted_state(double* dotz, 
                            double t,
                            const double* z, 
                            const params* parms)
{
  const double theta1 = z[0];
  const double theta2 = z[1];
  const double theta1dot = z[2];
  const double theta2dot = z[3];

  const double sn1 = std::sin(theta1);
  const double cs1 = std::cos(theta1);
  const double sn2 = std::sin(theta2);
  const double cs2 = std::cos(theta2);

  const double m1 = parms->m1;
  const double I1 = parms->I1;
  const double hL1 = parms->L1 / 2.0;

  const double m2 = parms->m2;
  const double I2 = parms->I2;
  const double hL2 = parms->L2 / 2.0;

  const double muA = parms->muA;
  const double muB = parms->muB;
  const double u = parms->u;
  const double g = parms->g;

  const double hL1sn1 = hL1 * sn1;
  const double hL1cs1 = hL1 * cs1;

  const double a11 = 1.0 / m1 + (hL1sn1 * hL1sn1) / I1;
  const double a12 = - (hL1sn1 * hL1cs1) / I1;
  const double a13 = 1.0 / m1 - (hL1sn1 * hL1sn1) / I1;
  const double a14 = -1.0 * a12;

  const double a21 = a12;
  const double a22 = 1.0 / m1 + (hL1cs1 * hL1cs1) / I1;
  const double a23 = a14;
  const double a24 = 1.0 / m1 - (hL1cs1 * hL1cs1) / I1;

  const double hL2sn2 = hL2 * sn2;
  const double hL2cs2 = hL2 * cs2;

  const double a31 = a13;
  const double a32 = a23;
  const double a33 = a11 + 1.0 / m2 + (hL2sn2 * hL2sn2) / I2;
  const double a34 = a12 - (hL2sn2 * hL2cs2) / I2;

  const double a41 = a14;
  const double a42 = a24;
  const double a43 = a34;
  const double a44 = a22 + 1.0 / m2 + (hL2cs2 * hL2cs2) / I2;

  const double mu1 = -muA * theta1dot - muB * (theta1dot - theta2dot);
  const double mu2 = muB * (theta1dot - theta2dot);

  const double sq_theta1dot = theta1dot * theta1dot;
  const double sq_theta2dot = theta2dot * theta2dot;

  const double bf1 = -hL1 * sq_theta1dot * cs1 - hL1 * sn1 * (u + mu1) / I1;
  const double bf2 = -hL1 * sq_theta1dot * sn1 + hL1 * cs1 * (u + mu1) / I1 + g;
  const double bf3 = hL1 * sq_theta1dot * cs1 + hL1 * sn1 * (u + mu1) / I1 + hL2 * sq_theta2dot * cs2 + hL2 * sn2 * (-u + mu2) / I2;
  const double bf4 = hL1 * sq_theta1dot * sn1 - hL1 * cs1 * (u + mu1) / I1 + hL2 * sq_theta2dot * sn2 - hL2 * cs2 * (-u + mu2) / I2;

  // AF = [a11,..,a14; a21,..,a24; a31,..,a34; a41,..,a44] is a symmetric matrix, with pos. diagonal
  const double AF[16] = {a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34, a41, a42, a43, a44};
  const double BF[4]  = {bf1, bf2, bf3, bf4};

  double Fab[4] = {0.0, 0.0, 0.0, 0.0}; // stores f = [fax;fay;fbx;fby];
  pivotfree_solve_4x4(AF, BF, Fab);

  const double theta1ddot = hL1 * (Fab[0] * sn1 - Fab[1] * cs1 - Fab[2] * sn1 + Fab[3] * cs1) / I1 + (u + mu1) / I1;
  const double theta2ddot = hL2 * (-Fab[2] * sn2 + Fab[3] * cs2) / I2 + (-u + mu2) / I2;

  dotz[0] = theta1dot;
  dotz[1] = theta2dot;
  dotz[2] = theta1ddot;
  dotz[3] = theta2ddot;
}

void calculate_CM_state(double* xy1, // [x1, y1, x1dot, y1dot]
                        double* xy2, // [x2, y2, x2dot, y2dot]
                        const double* z, 
                        const params* parms)
{
  const double theta1 = z[0];
  const double theta2 = z[1];
  const double theta1dot = z[2];
  const double theta2dot = z[3];

  const double sn1 = std::sin(theta1);
  const double cs1 = std::cos(theta1);
  const double sn2 = std::sin(theta2);
  const double cs2 = std::cos(theta2);

  const double hL1 = parms->L1 / 2.0;
  const double hL2 = parms->L2 / 2.0;

  const double x1 = cs1 * hL1;
  const double y1 = sn1 * hL1;

  const double x2 = 2.0 * x1 + cs2 * hL2;
  const double y2 = 2.0 * y1 + sn2 * hL2;

  const double x1dot = -1.0 * sn1 * hL1 * theta1dot;
  const double y1dot = 1.0 * cs1 * hL1 * theta1dot;

  const double x2dot = 2.0 * x1dot - 1.0 * sn2 * hL2 * theta2dot;
  const double y2dot = 2.0 * y1dot + 1.0 * cs2 * hL2 * theta2dot;

  xy1[0] = x1;
  xy1[1] = y1;
  xy1[2] = x1dot;
  xy1[3] = y1dot;

  xy2[0] = x2;
  xy2[1] = y2;
  xy2[2] = x2dot;
  xy2[3] = y2dot;
}

double total_mechanical_energy(const double* z,
                               const double* xy1,
                               const double* xy2,
                               const params* parms)
{
  const double Krot = 0.5 * parms->I1 * z[2] * z[2] + 
                      0.5 * parms->I2 * z[3] * z[3];

  const double Klin = 0.5 * parms->m1 * (xy1[2] * xy1[2] + xy1[3] * xy1[3]) + 
                      0.5 * parms->m2 * (xy2[2] * xy2[2] + xy2[3] * xy2[3]);

  const double Vgra = parms->m1 * parms->g * xy1[1] + 
                      parms->m2 * parms->g * xy2[1];
                      
  return (Krot + Klin + Vgra);
}

} // end namespace acrobot

/*
#if defined(HAVE_OCTAVE) || defined(MATLAB_MEX_FILE)
#endif
*/
