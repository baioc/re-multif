% Calculates the molecule concentration gradient of the network (using the
% Stochastic model) for given molecule concentrations R at instant t. I(t)
% should evaluate to the input signal at time t; Ha(X,k) and Hr(X,k) must
% compute the Hill functions for activation and repression of protein X given
% half-saturation constant k; and k_tl beta, P_tc, delta_x and delta_m are
% reaction constants extracted from the network. k_a and k_r are vectors with
% the binding affinity for each repressor protein + input. N is used as a matrix
% where each row represents noise at time t, with 12 columns for the coefficients.
function dR = stochastic(R, t, I, N, Ha, Hr, k_a, k_r, k_tl, beta, P_tc, delta_x, delta_m)

  % these are just coefficient indexes for the noise matrix
  a1_idx = 1;
  a2_idx = 2;
  a3_idx = 3;
  b1_idx = 4;
  b2_idx = 5;
  b3_idx = 6;
  c1_idx = 7;
  c2_idx = 8;
  c3_idx = 9;
  d1_idx = 10;
  d2_idx = 11;
  d3_idx = 12;

  % common subexpressions for coefficient formulas
  alpha = k_tl * beta / delta_m;
  gamma = k_tl * P_tc / delta_m;

  % coefficients
  a1 = alpha * Ha(I(t), k_a(5)) * Hr(R(2), k_r(2));
  a2 = gamma * Hr(R(3), k_r(3));
  a3 = delta_x * R(1);

  b1 = alpha * Ha(I(t), k_a(5)) * Hr(R(4), k_r(4));
  b2 = gamma * Hr(R(3), k_r(3)) * Hr(R(4), k_r(4));
  b3 = delta_x * R(2);

  c1 = alpha * Ha(I(t), k_a(5)) * Hr(R(4), k_r(4));
  c2 = gamma * Hr(R(1), k_r(1));
  c3 = delta_x * R(3);

  d1 = alpha * Ha(I(t), k_a(5)) * Hr(R(2), k_r(2));
  d2 = gamma * Hr(R(1), k_r(1)) * Hr(R(2), k_r(2));
  d3 = delta_x * R(4);

  V = 1e9; % "volumetric conversion factor"

  dR = zeros(size(R));

  dR(1) = (sqrt(a1)*N(t,a1_idx) + sqrt(a2)*N(t,a2_idx) - sqrt(a3)*N(t,a3_idx)) ...
          / sqrt(V) + a1 + a2 - a3;

  dR(2) = (sqrt(b1)*N(t,b1_idx) + sqrt(b2)*N(t,b2_idx) - sqrt(b3)*N(t,b3_idx)) ...
          / sqrt(V) + b1 + b2 - b3;

  dR(3) = (sqrt(c1)*N(t,c1_idx) + sqrt(c2)*N(t,c2_idx) - sqrt(c3)*N(t,c3_idx)) ...
          / sqrt(V) + c1 + c2 - c3;

  dR(4) = (sqrt(d1)*N(t,d1_idx) + sqrt(d2)*N(t,d2_idx) - sqrt(d3)*N(t,d3_idx)) ...
          / sqrt(V) + d1 + d2 - d3;

endfunction
