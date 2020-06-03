% Approximates the fractional saturation of a molecule X in the presence of
% negative-cooperative ligands with half-saturation k at n binding sites.
function h = hill_repression(X, n, k)

  h = 1 / (1 + (X.^n / k^n));

endfunction
