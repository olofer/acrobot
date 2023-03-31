% CLI usage: octave --no-window-system --eval "acrobot_rebuild_test"

isOctave = (exist('OCTAVE_VERSION', 'builtin') ~= 0);
assert(isOctave, 'script only supports Octave at the moment');
mexfilename = sprintf('acrobot_odefun.%s', mexext());

disp(sprintf('cleaning \"%s\"..', mexfilename));
autoload('acrobot_odefun', file_in_loadpath(mexfilename), "remove"); 
delete(mexfilename);

disp('building..');
mkoctfile --verbose --mex --strip acrobot_odefun.cpp;

try 
  pkg('load', 'control');
  hasControlPackage = true;
catch
  hasControlPackage = false;
end_try_catch

% Run the Octave simulator program & then resimulate using C++/MEX version
disp('testing (reference)..');
muA = 0.0125;
muB = 0.0250;
rep = acrobot_test(muA, muB);
P = rep.params;

fprintf(1, 'checking stationary points (has control package: %i)..\n', hasControlPackage);
Zeqs = [-pi/2, -pi/2, 0, 0;
        -pi/2,  pi/2, 0, 0; 
         pi/2, -pi/2, 0, 0;
         pi/2,  pi/2, 0, 0];
epfd = 1.0e-4;
dotZeqs = NaN(4, 4);
for i = 1:4
  dotZeqs(i, :) = acrobot_odefun(0, Zeqs(i, :), P, 0);
  assert(norm(dotZeqs(i, :)) < 1.0e-14, 'not stationary point');
  [Ai, Bi] = acrobot_odefun(epfd, Zeqs(i, :), P, 0);
  disp(sum(real(eig(Ai)) > 0)); % output should be: 0, 1, 1, 2
  if hasControlPackage
    % Confirm that each stationary point is controllable with state feedback
    assert(isctrb(Ai, Bi), 'not controllable stationary point');
  end
end

% Check equivalance to above simulation
Tref = rep.t;
Zref = rep.z(:, [3, 9, 6, 12]); % non-redundant state vector [theta1, theta2, theta1dot, theta2dot]
Uref = rep.u;

disp('testing (C++/MEX)..');
cpp_odefunc = @(t, a)(acrobot_odefun(t, a, P, interp1(Tref, Uref, t)));
ode_struct = odeset('RelTol', 1e-8, 'AbsTol', 1e-9);
[Tcpp, Zcpp] = ode45(cpp_odefunc, [0.0, 8.0], Zref(1, :), ode_struct);

err = Zref - interp1(Tcpp, Zcpp, Tref, 'spline');
relerr = max(abs(err)) ./ max(abs(Zref));
disp(relerr);

if max(relerr) < 1e-4
  disp('code seems to be working OK');
else
  warning('code may not be working OK');
end

if have_window_system() % Visual confirmation, if possible
  figure;
  plot(Tref, Zref, 'LineWidth', 3);
  hold on;
  plot(Tcpp, Zcpp, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
  xlabel('time');
  legend('theta1', 'theta2', 'dot(theta1)', 'dot(theta2)', '(C++/MEX)');
  grid on;
  title('Checking C++/MEX calculation');
end
