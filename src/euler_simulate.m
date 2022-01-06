% Simulates the model defined by function f, for a number of steps -- defined
% based on the given interval -- of dt seconds each, over an initial state X.
% f(t, X) should return the gradient Y at time t, where X and Y are vectors.
%
% Returns a vector T and a matrix Y, corresponding to each time step and protein
% concentrations at each of those steps.
function [T, Y] = euler_simulate(f, X, interval, dt)

  assert(size(X, 1) == 1);
  assert(interval(1) >= 0);
  assert(interval(2) > interval(1));
  assert(dt >= 1);

  m = ceil((interval(2) - interval(1)) / dt);
  n = length(X);
  Y = zeros(m, n);
  T = zeros(m, 1);

  for i = 1 : m
    Y(i,:) = X;
    t = (i-1) * dt;
    T(i,1) = interval(1) + t;
    dX = f(t, X);
    X += dX * dt;
  endfor

endfunction
