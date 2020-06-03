% Calculates the protein concentrations gradient of the network (using the QSSA
% model) for given repressor concentrations R at instant t. I(t) should evaluate
% to the input signal at time t; Ha(X,k) and Hr(X,k) must compute the Hill
% functions for activation and repression of protein X given half-saturation
% constant k; and alpha, gamma and delta are reaction constants extracted from
% the network. k_a and k_r are vectors with the binding affinity for each
% repressor protein + input.
function dR = qssa(R, t, I, Ha, Hr, k_a, k_r, alpha, gamma, delta)

  dR = zeros(size(R));

  dR(1) = alpha * Ha(I(t), k_a(5)) * Hr(R(2), k_r(2)) ...
          + gamma * Hr(R(3), k_r(3)) ...
          - delta * R(1);

  dR(2) = alpha * Ha(I(t), k_a(5)) * Hr(R(4), k_r(4)) ...
          + gamma * Hr(R(3), k_r(3)) * Hr(R(4), k_r(4)) ...
          - delta * R(2);

  dR(3) = alpha * Ha(I(t), k_a(5)) * Hr(R(4), k_r(4)) ...
          + gamma * Hr(R(1), k_r(1)) ...
          - delta * R(3);

  dR(4) = alpha * Ha(I(t), k_a(5)) * Hr(R(2), k_r(2)) ...
          + gamma * Hr(R(1), k_r(1)) * Hr(R(2), k_r(2)) ...
          - delta * R(4);

endfunction
