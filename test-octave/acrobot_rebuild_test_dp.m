isOctave = (exist('OCTAVE_VERSION', 'builtin') ~= 0);
assert(isOctave, 'script only supports Octave at the moment');
mexfilename = sprintf('acrobot_dpsolve.%s', mexext());

disp(sprintf('cleaning \"%s\"..', mexfilename));
autoload('acrobot_dpsolve', file_in_loadpath(mexfilename), "remove"); 
delete(mexfilename);

disp('building..');
mkoctfile --verbose --mex --strip acrobot_dpsolve.cpp;

P = struct;
P.L1 = 1.50;
P.M1 = 1.00;
P.L2 = 1.00;
P.M2 = 0.50;
P.I1 = P.M1 * (P.L1^2) / 12;
P.I2 = P.M2 * (P.L2^2) / 12;
P.g = 9.82;
P.muA = 0.0;
P.muB = 0.0;

disp('calling DP solver..');
%npts = [36, 34, 33, 35];
%npts = [25, 25, 50, 50];
npts = [54, 54, 65, 65];
itrs = 750;
deltat = 1.0e-2;

%acrobot_dpsolve(P, npts, itrs, deltat);

[V, A] = acrobot_dpsolve(P, npts, itrs, deltat);
%disp(size(V));
%disp(size(A));
%disp(prod(size(V)));

sum(V(:))
sum(A(:))

disp('dumping 4d tables to disk..');
Vf = single(V);
Af = single(A);
%save latest-dpsolve.mat npts itrs deltat P Vf Af;
save('-mat7-binary', 'latest-dpsolve.mat', 'npts', 'itrs', 'deltat', 'Vf', 'Af');
