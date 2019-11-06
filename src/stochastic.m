% ============================== OCTAVE SETUP ==================================
    clear;
    format long;
    pkg load statistics;
% ========================== MODEL REACTION RATES ==============================

    % translation rates (1/s)
    Kt_R1 = 6e-4;
    Kt_R2 = 6e-4;
    Kt_R3 = 6e-4;
    Kt_R4 = 6e-4;

    % mRNA degradation rates (1/s)
    Kd_mR1P1 = 2.5e-3;
    Kd_mR4P1 = 2.5e-3;
    Kd_mR2P2 = 2.5e-3;
    Kd_mR3P2 = 2.5e-3;
    Kd_mR1P3 = 2.5e-3;
    Kd_mR2P4 = 2.5e-3;
    Kd_mR3P5 = 2.5e-3;
    Kd_mR4P6 = 2.5e-3;

    % protein degradation rates (1/s)
    Kd_R1 = 4e-4;
    Kd_R2 = 4e-4;
    Kd_R3 = 4e-4;
    Kd_R4 = 4e-4;

    % maximum/unrepressed transcription rates (M/s)
    Kb_P1 = 4e-10;
    Kb_P2 = 4e-10;
    Kb_P3 = 4e-10;
    Kb_P4 = 4e-10;
    Kb_P5 = 4e-10;
    Kb_P6 = 4e-10;

    % Hill coefficients (scalar)
    Na_P1 = 1.3;
    Na_P2 = 1.3;
    Nr_R2P1 = 1.3;
    Nr_R4P2 = 1.3;
    Nr_R3P3 = 1.3;
    Nr_R3P4 = 1.3;
    Nr_R4P4 = 1.3;
    Nr_R1P5 = 1.3;
    Nr_R1P6 = 1.3;
    Nr_R2P6 = 1.3;

    % Hill activation constants (M)
    Ka_P1 = 2e-8;
    Ka_P2 = 2e-8;

    % Hill repression constants (M)
    Kr_R2P1 = 6e-10;
    Kr_R4P2 = 6e-10;
    Kr_R3P3 = 6e-10;
    Kr_R3P4 = 6e-10;
    Kr_R4P4 = 6e-10;
    Kr_R1P5 = 6e-10;
    Kr_R1P6 = 6e-10;
    Kr_R2P6 = 6e-10;


% ============================ INPUT PARAMETERS ================================

    % time settings
    dt = 60;                     % 1 minute timesteps
    simulation = 0 : dt : 4.0e5; % simulation length control

    % input vector
    Is = (6 .- 25*cos(simulation .* 2*pi/0.9e5) .+ 25) * 1e-9; % sinusoidal

    % initial state concentrations (M)
    R1 = 50e-9;
    R2 = 50e-9;
    R3 = 0e-9;
    R4 = 0e-9;


% =============================== SIMULATION ===================================

    % preallocate arrays to store concentrations over time
    R1s = zeros(1, length(simulation));
    R2s = zeros(1, length(simulation));
    R3s = zeros(1, length(simulation));
    R4s = zeros(1, length(simulation));

    % precompute loop invariants
    Ka_P1 ^= Na_P1;
    Ka_P2 ^= Na_P2;
    Kr_R2P1 ^= Nr_R2P1;
    Kr_R4P2 ^= Nr_R4P2;
    Kr_R3P3 ^= Nr_R3P3;
    Kr_R3P4 ^= Nr_R3P4;
    Kr_R4P4 ^= Nr_R4P4;
    Kr_R1P5 ^= Nr_R1P5;
    Kr_R1P6 ^= Nr_R1P6;
    Kr_R2P6 ^= Nr_R2P6;

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
        AP1 = Ha(Is(t), Na_P1, Ka_P1);
        AP2 = Ha(Is(t), Na_P2, Ka_P2);
        RP1 = Hr(R2, Nr_R2P1, Kr_R2P1);
        RP2 = Hr(R4, Nr_R4P2, Kr_R4P2);

        % simplified coefficients
        a1 = (Kt_R1 * Kb_P1 / Kd_mR1P1) * AP1 * RP1;
        a2 = (Kt_R1 * Kb_P3 / Kd_mR1P3) * Hr(R3, Nr_R3P3, Kr_R3P3);
        a3 = Kd_R1 * R1;
        b1 = (Kt_R2 * Kb_P2 / Kd_mR2P2) * AP2 * RP2;
        b2 = (Kt_R2 * Kb_P4 / Kd_mR2P4) * Hr(R3, Nr_R3P4, Kr_R3P4) * Hr(R4, Nr_R4P4, Kr_R4P4);
        b3 = Kd_R2 * R2;
        c1 = (Kt_R3 * Kb_P2 / Kd_mR3P2) * AP2 * RP2;
        c2 = (Kt_R3 * Kb_P5 / Kd_mR3P5) * Hr(R1, Nr_R1P5, Kr_R1P5);
        c3 = Kd_R3 * R3;
        d1 = (Kt_R4 * Kb_P1 / Kd_mR4P1) * AP1 * RP1;
        d2 = (Kt_R4 * Kb_P6 / Kd_mR4P6) * Hr(R1, Nr_R1P6, Kr_R1P6) * Hr(R2, Nr_R2P6, Kr_R2P6);
        d3 = Kd_R4*R4;

        V = 1e9;

        noise = stdnormal_rnd(1)

        % compute variation
        dR1dt = a1 + a2 - a3 + ((1/sqrt(V)) * ((sqrt(a1) * noise) + (sqrt(a2) *
        noise) - (sqrt(a3) * noise)));
        dR2dt = b1 + b2 - b3 + ((1/sqrt(V)) * ((sqrt(b1) * noise) + (sqrt(b2) *
        noise) - (sqrt(b3) * noise)));
        dR3dt = c1 + c2 - c3 + ((1/sqrt(V)) * ((sqrt(c1) * noise) + (sqrt(c2) *
        noise) - (sqrt(c3) * noise)));
        dR4dt = d1 + d2 - d3 + ((1/sqrt(V)) * ((sqrt(d1) * noise) + (sqrt(d2) *
        noise) - (sqrt(d3) * noise)));

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
    simulation /= timescale;
    Is  /= quantscale;
    R1s /= quantscale;
    R2s /= quantscale;
    R3s /= quantscale;
    R4s /= quantscale;

    % model suplot
    subplot(2, 1, 1);
    plot(simulation,R1s,'-m;R1;', simulation,R2s,'-k;R2;', simulation,R3s,'-r;R3;', simulation,R4s,'-g;R4;');
    xlabel("Time (10^5 seconds)");
    ylabel("Concentration (nM)");
    title("Reduced Model");

    % input subplot
    subplot(2, 1, 2);
    plot(simulation, Is, 'b;I;');
    xlabel("Time (10^5 seconds)");
    ylabel("Concentration (nM)");

    hold off;
    a = input("\nPress enter to exit ");
