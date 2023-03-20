function rep = acrobot_test(muA, muB)
%
% function rep = acrobot_test(muA, muB)
%
% Test program to verify basic governing equations for the Acrobot.
%
% Integrate ODE in redundant configuration space to confirm that reaction forces 
% are correctly solved for such that they ensure the joint constraints are
% maintained (up to some reasonable accuracy). Also confirm that energy is conserved.
% Confirmation (for now) is simply done visually by checking the plots generated.
%
% The friction parameters muA, muB > 0 are optional (they default to zero).
%

P = struct;

P.L1 = 1.50;
P.M1 = 1.00;
P.L2 = 1.00;
P.M2 = 0.50;

% thin rod inertia w.r.t CM
P.I1 = P.M1 * (P.L1^2) / 12;
P.I2 = P.M2 * (P.L2^2) / 12;

P.g = 9.82;

if nargin < 2, muB = 0; end 
if nargin < 1, muA = 0; end
assert(muA >= 0 && muB >= 0, 'friction coefficients must be non-negative');
P.muA = muA;
P.muB = muB;

S = struct;

% Initial configuration
S.theta1 = 0.0;
S.theta1dot = 0.0;
S.theta2 = 0.0;
S.theta2dot = 0.0;

S.u = 0.0; % applied joint torque

zinit = make_z_vector(S, P);

odefunc = @(t, a)(redundant_ode_system(t, a, P));

ode_struct = odeset('RelTol', 1e-8, 'AbsTol', 1e-9);
[t, z] = ode45(odefunc, [0.0, 8.0], zinit, ode_struct);

[Kin, Pot] = calc_energies(z, P);

rep = struct;
rep.creator = mfilename();

rep.t = t;
rep.z = z;

rep.K = Kin;
rep.V = Pot;

rep.u = NaN(size(rep.t));
for k = 1:numel(rep.t)
  rep.u(k) = make_torque_input(rep.t(k));
end

rep.Pu = z(:, 6) .* rep.u - z(:, 12) .* rep.u;
rep.params = P;

if nargout == 0
  figure;
  hold on;
  plot(z(:, 1), z(:, 2), 'b-');
  plot(z(:, 7), z(:, 8), 'r-');

  % plot final arm position
  zf = z(end, :);
  p0 = [zf(1) - cos(zf(3)) * P.L1 / 2, zf(2) - sin(zf(3)) * P.L1 / 2];
  p1 = [zf(1) + cos(zf(3)) * P.L1 / 2, zf(2) + sin(zf(3)) * P.L1 / 2];
  line([p0(1), p1(1)], [p0(2), p1(2)], 'Color', [0.5, 0, 1], 'LineWidth', 2);

  p2 = [zf(7) - cos(zf(9)) * P.L2 / 2, zf(8) - sin(zf(9)) * P.L2 / 2];
  p3 = [zf(7) + cos(zf(9)) * P.L2 / 2, zf(8) + sin(zf(9)) * P.L2 / 2];
  line([p2(1), p3(1)], [p2(2), p3(2)], 'Color', [1, 0, 0.5], 'LineWidth', 2);

  legend('Arm 1', 'Arm 2');
  grid on;
  axis equal;
  xlabel('x');
  ylabel('y');
  title('Center-of-mass movements');

  figure;
  plot(rep.t, [rep.K, rep.V, rep.K + rep.V, cumtrapz(rep.t, rep.Pu)], 'LineWidth', 2);
  legend('kinetic', 'potential', 'total', 'input');
  grid on;
  xlabel('time');
  ylabel('Energy [J]');
  title('Conservation of mechanical energy');
end

end

function [FA, FB, ddot1, ddot2] = calc_reactions(P, S)
  A11 = [1, 0, sin(S.theta1) * P.L1 / 2; 
         0, 1, -1 * cos(S.theta1) * P.L1 / 2];
  b1 = (-1 * S.theta1dot^2 * P.L1 / 2) * [cos(S.theta1); sin(S.theta1)];

  A21 = [1, 0, -1 * sin(S.theta1) * P.L1 / 2;
         0, 1, cos(S.theta1) * P.L1 / 2];
  A22 = [-1, 0, -1 * sin(S.theta2) * P.L2 / 2;
         0, -1, cos(S.theta2) * P.L2 / 2];
  b2 = [S.theta1dot^2 * cos(S.theta1) * P.L1 / 2 + S.theta2dot^2 * cos(S.theta2) * P.L2 / 2;
        S.theta1dot^2 * sin(S.theta1) * P.L1 / 2 + S.theta2dot^2 * sin(S.theta2) * P.L2 / 2];

  M1 = diag([P.M1, P.M1, P.I1]);
  B1 = [1, 0, 1, 0;
        0, 1, 0, 1;
        sin(S.theta1) * P.L1 / 2, -1 * cos(S.theta1) * P.L1 / 2, -1 * sin(S.theta1) * P.L1 / 2, cos(S.theta1) * P.L1 / 2];
  u1 = [0; -P.M1 * P.g; S.u - P.muA * S.theta1dot - P.muB * (S.theta1dot - S.theta2dot)];

  M2 = diag([P.M2, P.M2, P.I2]);
  B2 = [0, 0, -1, 0;
        0, 0, 0, -1;
        0, 0, -1 * sin(S.theta2) * P.L2 / 2, cos(S.theta2) * P.L2 / 2];
  u2 = [0; -P.M2 * P.g; -S.u + P.muB * (S.theta1dot - S.theta2dot)];

  % Assemble final equation for reaction forces & solve it. 

  A = [A11 * (M1\B1); 
       A21 * (M1\B1) + A22 * (M2\B2)];
  b = [b1 - A11 * (M1\u1);
       b2 - A21 * (M1\u1) - A22 * (M2\u2)];
  F = A\b;

  % Return results including all net accelerations.

  FA = F(1:2);
  FB = F(3:4);
  ddot1 = M1 \ (B1 * F + u1);
  ddot2 = M2 \ (B2 * F + u2);
end

% Define full vector z from [theta1, theta1dot, theta2, theta2dot]
function z = make_z_vector(S, P)

  x1 = 0 + cos(S.theta1) * P.L1 / 2;
  y1 = 0 + sin(S.theta1) * P.L1 / 2;
  x2 = 0 + cos(S.theta1) * P.L1 + cos(S.theta2) * P.L2 / 2;
  y2 = 0 + sin(S.theta1) * P.L1 + sin(S.theta2) * P.L2 / 2;

  x1dot = -1 * (P.L1 / 2) * sin(S.theta1) * S.theta1dot;
  y1dot = (P.L1 / 2) * cos(S.theta1) * S.theta1dot;
  x2dot = -1 * P.L1* sin(S.theta1) * S.theta1dot - (P.L2 / 2) * sin(S.theta2) * S.theta2dot;
  y2dot = P.L1 * cos(S.theta1) * S.theta1dot + (P.L2 / 2) * cos(S.theta2) * S.theta2dot;

  z1 = [x1; y1; S.theta1; x1dot; y1dot; S.theta1dot];
  z2 = [x2; y2; S.theta2; x2dot; y2dot; S.theta2dot];
  z = [z1; z2];
end

function dz = redundant_ode_system(t, z, P)
  assert(numel(z) == 12); % z = [x1, y1, theta1, x1dot, y1dot, theta1dot,
                          %      x2, y2, theta2, x2dot, y2dot, theta2dot]
  state = struct;
  state.u = make_torque_input(t);
  state.theta1 = z(3);
  state.theta1dot = z(6);
  state.theta2 = z(9);
  state.theta2dot = z(12);

  [FA, FB, ddot1, ddot2] = calc_reactions(P, state);

  dz = [z(4:6); ddot1; z(10:12); ddot2];
end

function [Kin, Pot] = calc_energies(z, P)
  assert(size(z, 2) == 12);
  Klin = 0.5 * P.M1 * (z(:, 4).^2 + z(:, 5).^2) + 0.5 * P.M2 * (z(:, 10).^2 + z(:, 11).^2);
  Krot = 0.5 * P.I1 * z(:, 6).^2 + 0.5 * P.I2 * z(:, 12).^2;
  Kin = Klin + Krot;
  Pot = P.M1 * P.g * z(:, 2) + P.M2 * P.g * z(:, 8);
end

function u = make_torque_input(t)
  if t > 5.0 && t < 6.0
    u = 1.25;
  else 
    u = 0.00;
  end
end 
