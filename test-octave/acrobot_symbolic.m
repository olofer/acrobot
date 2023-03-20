%
% Script to symbolically derive coefficients in the system of equations
% for the reaction forces: AF * f = BF, where f = [Fax;Fay;Fbx;Fby].
%
% AF only depends on configuration:
%   q = (th1, th2);
%
% whereas BF depends on (q, qdot) and u.
%

pkg load symbolic;

th1 = sym('th1');
th1d = sym('th1d');

th2 = sym('th2');
th2d = sym('th2d');

L1 = sym('L1');
I1 = sym('I1');
m1 = sym('m1');

L2 = sym('L2');
I2 = sym('I2');
m2 = sym('m2');

g = sym('g');
u = sym('u');

muA = sym('muA');
muB = sym('muB');

A11 = [1, 0, (L1/2)*sin(th1);
       0, 1, (-L1/2)*cos(th1)];

b1 = [cos(th1); sin(th1)] * (-L1/2) * th1d^2;

A21 = [1, 0, (-L1/2)*sin(th1);
       0, 1, (L1/2)*cos(th1)];

A22 = [-1, 0, (-L2/2)*sin(th2);
       0, -1, (L2/2)*cos(th2)];

b2 = [(L1/2)*th1d^2*cos(th1) + (L2/2)*th2d^2*cos(th2);
      (L1/2)*th1d^2*sin(th1) + (L2/2)*th2d^2*sin(th2)];

M1 = diag([m1; m1; I1]);

B1 = [sym(1), 0, 1, 0;
      sym(0), 1, 0, 1;
      (L1/2)*sin(th1), (-L1/2)*cos(th1), (-L1/2)*sin(th1), (L1/2)*cos(th1)];

u1 = [0; -m1*g; u - muA * th1d - muB * (th1d - th2d)];

M2 = diag([m2; m2; I2]);

B2 = [sym(0), 0, -1, 0;
      sym(0), 0, 0, -1;
      0, 0, (-L2/2)*sin(th2), (L2/2)*cos(th2)];

u2 = [0; -m2*g; -u + muB * (th1d - th2d)];

AF = [A11*inv(M1)*B1; 
      A21*inv(M1)*B1 + A22*inv(M2)*B2];
BF = [b1 - A11*inv(M1)*u1;
      b2 - A21*inv(M1)*u1 - A22*inv(M2)*u2];

% Now the f = AF\BF; in principle (do not try).
%
% and then M1*ddot1 = B1*f + u1 
%          M2*ddot2 = B2*f + u2
%
% where ddot1 = dd(x1,y1,th1) and similarly for 2 in place of 1.
% with dd = double-dot. So only the last row of the above equations
% are needed to close the system for the (q,qdot) state variable.

disp(AF - AF.'); % should show a 4x4 matrix of zeros

% Evaluate all coefficients in the AF matrix; print out the diagonal elements

a11 = simplify(AF(1, 1));
a12 = simplify(AF(1, 2));
a13 = simplify(AF(1, 3));
a14 = simplify(AF(1, 4));

pretty(a11);

a21 = simplify(AF(2, 1));
a22 = simplify(AF(2, 2));
a23 = simplify(AF(2, 3));
a24 = simplify(AF(2, 4));

pretty(a22);

a31 = simplify(AF(3, 1));
a32 = simplify(AF(3, 2));
a33 = simplify(AF(3, 3));
a34 = simplify(AF(3, 4));

pretty(a33);

a41 = simplify(AF(4, 1));
a42 = simplify(AF(4, 2));
a43 = simplify(AF(4, 3));
a44 = simplify(AF(4, 4));

pretty(a44);

% Use these for numerical evaluation (if needed)

symbolOrderAF = symvar(AF);
matlabFunctionAF = matlabFunction(AF);

symbolOrderBF = symvar(BF);
matlabFunctionBF = matlabFunction(BF);

