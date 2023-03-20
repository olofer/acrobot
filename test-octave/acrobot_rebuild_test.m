isOctave = (exist('OCTAVE_VERSION', 'builtin') ~= 0);
assert(isOctave, 'script only supports Octave at the moment');
mexfilename = sprintf('acrobot_odefun.%s', mexext());

disp(sprintf('cleaning \"%s\"..', mexfilename));
autoload('acrobot_odefun', file_in_loadpath(mexfilename), "remove"); 
delete(mexfilename);

disp('building..');
mkoctfile --verbose --mex --strip acrobot_odefun.cpp;

% Run the Octave simulator program & then resimulate using C++/MEX version
disp('testing (reference)..');
muA = 0.0125;
muB = 0.0250;
rep = acrobot_test(muA, muB);
P = rep.params;

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
