% Evaluates a square wave signal at instant t. The signal has a period T,
% amplitude of A, an offset level DC and its duty cycle is duty.
function y = square_wave(t, T, A = 1, DC = 0, duty = 0.5)

  assert(T > 0);
  assert(0 <= duty && duty <= 1.0);

  if mod(t, T) > duty * T
    y = A + DC;
  else
    y = DC;
  endif

endfunction
