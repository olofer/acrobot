%
% Load "latest-dpsolve.mat" and run a simulation
%
% USAGE:
%   octave --eval "acrobot_check_dp" --persist
%

show_example_slice = false;
use_dp_feebdack = true;

DP = load('newest-dpsolve.mat');
disp(fieldnames(DP));

disp(size(DP.Af));
assert(length(DP.g1) == size(DP.Af, 1));
assert(length(DP.g2) == size(DP.Af, 2));
assert(length(DP.g3) == size(DP.Af, 3));
assert(length(DP.g4) == size(DP.Af, 4));

disp('ulevels:');
disp(unique(DP.Af(:)));

disp(fieldnames(DP.P));

if show_example_slice
  imid3 = round(size(DP.Vf, 3) / 2);
  imid4 = round(size(DP.Vf, 4) / 2);

  figure;
  imagesc(DP.g1, DP.g2, DP.Vf(:, :, imid3, imid4));
  axis xy;
  axis equal;
  xlabel('theta1');
  ylabel('theta2');
  title('Value function slice');

  figure;
  imagesc(DP.g1, DP.g2, DP.Af(:, :, imid3, imid4));
  axis xy;
  axis equal;
  xlabel('theta1');
  ylabel('theta2');
  title('Action function slice');
end

dtsim = 2e-3; %DP.deltat; %5e-3;
N = 1e3;
Z = NaN(N, 4);
U = NaN(N, 1);
T = dtsim * (0:(N - 1));
%Z(1, :) = [-pi/2 - 1e-3, -pi/2 + 1e-3, 0, 0];
%Z(1, :) = [pi/2 - 1e-3, pi/2 + 1e-3, 0, 0];
Z(1, :) = [0, 0, 0, 0];
%Z(1, :) = [0, 0, 2, 3];
fprintf(1, 'simulating (dt=%f) for %i steps (DP dt=%f) ... \n', dtsim, N, DP.deltat);
for t = 1:(N - 1)
  zt = Z(t, :);
  if use_dp_feebdack
    zi = zt;
    while zi(1) < 0
      zi(1) = zi(1) + 2 * pi;
    end
    while zi(1) > 2 * pi
      zi(1) = zi(1) - 2 * pi;
    end
    while zi(2) < 0
      zi(2) = zi(2) + 2 * pi;
    end
    while zi(2) > 2 * pi
      zi(2) = zi(2) - 2 * pi;
    end
    ut = interpn(DP.g1, DP.g2, DP.g3, DP.g4, DP.Af, zi(1), zi(2), zi(3), zi(4), 'nearest', 0.0);
    ut = double(ut);
  else
    ut = 0; 
  end
  U(t) = ut;
  ztdot = acrobot_odefun(T(t), zt, DP.P, ut);
  zpred = zt + dtsim * ztdot;
  zpreddot = acrobot_odefun(T(t + 1), zpred, DP.P, ut);
  znext = zt + 0.5 * dtsim * (ztdot + zpreddot);
  Z(t + 1, :) = znext;
end

figure;
stairs(T, [Z, U], 'LineWidth', 2);
legend('th1 [rad]', 'th2 [rad]', 'th1d [rad/s]', 'th2d [rad/s]', 'u [Nm]');
xlabel('time [sec]');
hold on;
for ii = -1:1
  line([T(1), T(end)], 2 * pi * ii + [pi/2, pi/2], 'Color', 'c', 'LineStyle', '--');
  line([T(1), T(end)], 2 * pi * ii + [-pi/2, -pi/2], 'Color', 'y', 'LineStyle', '-.');
end 
line([T(1), T(end)], [0, 0], 'Color', 'k', 'LineStyle', ':');
