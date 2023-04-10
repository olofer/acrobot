isOctave = (exist('OCTAVE_VERSION', 'builtin') ~= 0);
assert(isOctave, 'script only supports Octave at the moment');
mexfilename = sprintf('acrobot_dpsolve.%s', mexext());

disp(sprintf('cleaning \"%s\"..', mexfilename));
autoload('acrobot_dpsolve', file_in_loadpath(mexfilename), 'remove'); 
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
npts = [75, 75, 80, 80]; 
itrs = 1200;
deltat = 8e-3;
tdisc = 5.0;

%acrobot_dpsolve(P, npts, itrs, deltat, tdisc);
%return;

[V, A, g1, g2, g3, g4] = acrobot_dpsolve(P, npts, itrs, deltat, tdisc);

disp(unique(A(:)));
disp([sum(V(:)), min(V(:)), max(V(:)), mean(V(:))]);

disp('dumping 4d tables to disk..');
Vf = single(V);
Af = single(A);
save('-mat7-binary', ...
     'newest-dpsolve.mat', ...
     'npts', 'itrs', 'deltat', 'Vf', 'Af', 'P', 'g1', 'g2', 'g3', 'g4');
