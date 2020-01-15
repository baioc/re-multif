% ============================== OCTAVE SETUP ==================================
    clear;
    format long;


% ========================== MODEL REACTION RATES ==============================

    Kt = 6e-4;     % translation rate (1/s)
    Kd_m = 2.5e-3; % mRNA degradation rate (1/s)
    Kd_p = 4e-4;   % protein degradation rate (1/s)
    N = 1.3;       % Hill coefficient (scalar)
    Ka = 2e-8;     % Hill activation constant (M)
    Kr = 6e-10;    % Hill repression constant (M)
    Kb = 4e-10;    % maximum/unrepressed transcription rate (M/s)


% ============================ INPUT PARAMETERS ================================

    % time settings
    dt = 60;                     % 1 minute timesteps
    simulation = 0 : dt : 7.0e5; % simulation length control

    % input vector
    period = 1.6e5;
    amplitude = 50e-9;
    minimum = 0e-9;
    duty_cycle = 0.5;
    Is = zeros(1, length(simulation));
    % square wave
    for t = 1 : length(Is)
        if mod(t*dt, period) > duty_cycle*period
            Is(t) = amplitude;
        else
            Is(t) = minimum;
        endif
    endfor

    % initial state concentrations (M)
    R1 = R2 = 50e-9;
    R3 = R4 = 0e-9;


% =============================== SIMULATION ===================================

    % preallocate arrays to store concentrations over time
    R1s = zeros(1, length(simulation));
    R2s = zeros(1, length(simulation));
    R3s = zeros(1, length(simulation));
    R4s = zeros(1, length(simulation));

    % precompute loop invariants and reduced coefficients
    transcription = Kt * Kb / Kd_m;
    Kan = Ka ^ N;
    Krn = Kr ^ N;

    % abstractions for activation & repression Hill-function
    % these consider K^n is being passed in already precomputed
    function y = Ha(X, n, Kn)
        Xn = X^n;
        y = Xn / (Kn + Xn);
    endfunction
    function y = Hr(X, n, Kn)
        y = 1 / (1 + (X^n / Kn));
    endfunction

    % actual simulation
    for t = 1 : length(simulation)

        % store current state
        R1s(t) = R1;
        R2s(t) = R2;
        R3s(t) = R3;
        R4s(t) = R4;

        % common subexpression optimization (used below)
        HaI = Ha(Is(t), N, Kan);
        HrR1 = Hr(R1, N, Krn);
        HrR2 = Hr(R2, N, Krn);
        HrR3 = Hr(R3, N, Krn);
        HrR4 = Hr(R4, N, Krn);

        % compute variation
        dR1dt = transcription * (HaI*HrR2 + HrR3) - Kd_p*R1;
        dR2dt = transcription * (HaI*HrR4 + HrR3*HrR4) - Kd_p*R2;
        dR3dt = transcription * (HaI*HrR4 + HrR1) - Kd_p*R3;
        dR4dt = transcription * (HaI*HrR2 + HrR1*HrR2) - Kd_p*R4;

        % apply state changes
        R1 += dR1dt * dt;
        R2 += dR2dt * dt;
        R3 += dR3dt * dt;
        R4 += dR4dt * dt;

    endfor


% ================================ RESULTS =====================================

    figure;
    hold on;

    % scale data for easier visualization
    timescale = 1e5;
    quantscale = 1e-9;
    x = simulation / timescale;
    yI = Is / quantscale;
    yR1 = R1s / quantscale;
    yR2 = R2s / quantscale;
    yR3 = R3s / quantscale;
    yR4 = R4s / quantscale;

    % model subplot
    subplot(2, 1, 1);
    plot(x,yR1,'--m;R1;', x,yR2,':k;R2;', x,yR3,'-r;R3;', x,yR4,'-.g;R4;');
    pbaspect([1 0.334 1]);
    xlabel("Time (10^5 seconds)");
    ylabel("Concentration (nM)");

    % input subplot
    subplot(2, 1, 2);
    plot(x, yI, 'b;I;');
    axis([-Inf, +Inf, 0, amplitude*1.23/quantscale]);
    pbaspect([1 0.334 1]);
    xlabel("Time (10^5 seconds)");
    ylabel("Concentration (nM)");

    hold off;
    print('freqdiv-square_.pdf'); % put in the folder the script is run from
