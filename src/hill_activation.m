% Approximates the fractional saturation of a molecule X in the presence of
% positive-cooperative ligands with half-saturation k at n binding sites.
function h = hill_activation(X, n, k)

  Xn = X .^ n;
  h = Xn / (k^n + Xn);

endfunction
