% Calculates the molecule concentration gradient of the network (using the Full
% model) for given molecule concentrations R at instant t. I(t) should evaluate
% to the input signal at time t; Ha(X,k) and Hr(X,k) must compute the Hill
% functions for activation and repression of protein X given half-saturation
% constant k; and k_tl beta, P_tc, delta_x and delta_m are reaction constants
% extracted from the network. k_a and k_r are vectors with the binding affinity
% for each repressor protein + input.
function dR = full_ode(R, t, I, Ha, Hr, k_a, k_r, k_tl, beta, P_tc, delta_x, delta_m)

  % these are just indexes for the mRNAs
  R1P1 = 1;
  R4P1 = 2;
  R2P2 = 3;
  R3P2 = 4;
  R1P3 = 5;
  R2P4 = 6;
  R3P5 = 7;
  R4P6 = 8;

  m = R(5:end); % mRNAs
  R = R(1:4);   % proteins

  dR = zeros(size(R));
  dm = zeros(size(m));

  % mRNA equations
  dm(R1P1) = beta * Ha(I(t), k_a(5)) * Hr(R(2), k_r(2)) - delta_m * m(R1P1);

  dm(R4P1) = beta * Ha(I(t), k_a(5)) * Hr(R(2), k_r(2)) - delta_m * m(R4P1);

  dm(R2P2) = beta * Ha(I(t), k_a(5)) * Hr(R(4), k_r(4)) - delta_m * m(R2P2);

  dm(R3P2) = beta * Ha(I(t), k_a(5)) * Hr(R(4), k_r(4)) - delta_m * m(R3P2);

  dm(R1P3) = P_tc * Hr(R(3), k_r(3))                    - delta_m * m(R1P3);

  dm(R2P4) = P_tc * Hr(R(3), k_r(3)) * Hr(R(4), k_r(4)) - delta_m * m(R2P4);

  dm(R3P5) = P_tc * Hr(R(1), k_r(1))                    - delta_m * m(R3P5);

  dm(R4P6) = P_tc * Hr(R(1), k_r(1)) * Hr(R(2), k_r(2)) - delta_m * m(R4P6);

  % protein equations
  dR(1) = k_tl * (m(R1P1) + m(R1P3)) - delta_x * R(1);
  dR(2) = k_tl * (m(R2P2) + m(R2P4)) - delta_x * R(2);
  dR(3) = k_tl * (m(R3P2) + m(R3P5)) - delta_x * R(3);
  dR(4) = k_tl * (m(R4P1) + m(R4P6)) - delta_x * R(4);

  % concatenate result
  dR = [dR dm];

endfunction
