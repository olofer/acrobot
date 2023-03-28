%
% Load "latest-dpsolve.mat" and run a simulation
%
% USAGE:
%   octave --eval "acrobot_check_dp" --persist
%

show_example_slice = false;

DP = load('latest-dpsolve.mat');
disp(fieldnames(DP));

disp(size(DP.Af));
assert(length(DP.g1) == size(DP.Af, 1));
assert(length(DP.g2) == size(DP.Af, 2));
assert(length(DP.g3) == size(DP.Af, 3));
assert(length(DP.g4) == size(DP.Af, 4));

disp('ulevels:');
disp(unique(DP.Af(:)));

if show_example_slice
  figure;
  imagesc(DP.g1, DP.g2, DP.Vf(:, :, 32, 34));
  axis xy;
  axis equal;
  xlabel('theta1');
  ylabel('theta2');
  title('Value function slice');

  figure;
  imagesc(DP.g1, DP.g2, DP.Af(:, :, 32, 34));
  axis xy;
  axis equal;
  xlabel('theta1');
  ylabel('theta2');
  title('Action function slice');
end

% TODO: run simulation !

