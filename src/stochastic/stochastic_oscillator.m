% ============================== OCTAVE SETUP ==================================
    clear;
    format long;


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

    % basal/leakage transcription rates (M/s)
    Kc_P1 = 0;
    Kc_P2 = 0;
    Kc_P3 = 0;
    Kc_P4 = 0;
    Kc_P5 = 0;
    Kc_P6 = 0;

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
    simulation = 0 : dt : 3.5e5; % simulation length control

    % input vector
    Is = zeros(1, length(simulation));

    % reaction initial conditions for each run
    multiI   = 1e-9 * 5; % constant input concentration
    multiR12 = 1e-9 * 50;
    multiR34 = 1e-9 * 0;


% ================================ MULTISIM ====================================

    % precompute loop invariants and reduced coefficients
    a1 = Kt_R1 * Kb_P1 / Kd_mR1P1;
    a2 = Kt_R1 * Kb_P3 / Kd_mR1P3;
    b1 = Kt_R2 * Kb_P2 / Kd_mR2P2;
    b2 = Kt_R2 * Kb_P4 / Kd_mR2P4;
    c1 = Kt_R3 * Kb_P2 / Kd_mR3P2;
    c2 = Kt_R3 * Kb_P5 / Kd_mR3P5;
    d1 = Kt_R4 * Kb_P1 / Kd_mR4P1;
    d2 = Kt_R4 * Kb_P6 / Kd_mR4P6;
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

    % for each set of parameters, make a simulation run and display its results

    % set initial state concentrations
    Is(:) = multiI;
    R1 = multiR12;
    R2 = multiR12;
    R3 = multiR34;
    R4 = multiR34;

    % preallocate arrays to store concentrations over time
    R1s = zeros(1, length(simulation));
    R2s = zeros(1, length(simulation));
    R3s = zeros(1, length(simulation));
    R4s = zeros(1, length(simulation));

    % generate Gaussian noise with zero mean and variance one
    random_seed = 73544911520192;
    randn('seed', random_seed);

    noise_a1 = randn(length(simulation), 1) / 50;
    noise_a2 = randn(length(simulation), 1) / 50;
    noise_a3 = randn(length(simulation), 1) / 50;
    noise_b1 = randn(length(simulation), 1) / 50;
    noise_b2 = randn(length(simulation), 1) / 50;
    noise_b3 = randn(length(simulation), 1) / 50;
    noise_c1 = randn(length(simulation), 1) / 50;
    noise_c2 = randn(length(simulation), 1) / 50;
    noise_c3 = randn(length(simulation), 1) / 50;
    noise_d1 = randn(length(simulation), 1) / 50;
    noise_d2 = randn(length(simulation), 1) / 50;
    noise_d3 = randn(length(simulation), 1) / 50;
    

    V = 1e9;
    isqrV = 1 / sqrt(V);

    % actual simulation
    for t = 1 : length(simulation)

        % store current state
        R1s(t) = R1;
        R2s(t) = R2;
        R3s(t) = R3;
        R4s(t) = R4;

        % common subexpression optimization (used below)
        activationP1 = Ha(Is(t), Na_P1, Ka_P1);
        repressionP1 = Hr(R2, Nr_R2P1, Kr_R2P1);
        activationP2 = Ha(Is(t), Na_P2, Ka_P2);
        repressionP2 = Hr(R4, Nr_R4P2, Kr_R4P2);

        a1_prime = a1 * activationP1 * repressionP1;
        a2_prime = a2 * Hr(R3, Nr_R3P3, Kr_R3P3);
        a3_prime = Kd_R1 * R1;
        b1_prime = b1 * activationP2 * repressionP2;
        b2_prime = b2 * Hr(R3, Nr_R3P4, Kr_R3P4) * Hr(R4, Nr_R4P4, Kr_R4P4);
        b3_prime = Kd_R2 * R2;
        c1_prime = c1 * activationP2 * repressionP2;
        c2_prime = c2 * Hr(R1, Nr_R1P5, Kr_R1P5);
        c3_prime = Kd_R3 * R3;
        d1_prime = d1 * activationP1 * repressionP1;
        d2_prime = d2 * Hr(R1, Nr_R1P6, Kr_R1P6) * Hr(R2, Nr_R2P6, Kr_R2P6);
        d3_prime = Kd_R4 * R4;

        % compute variation
        dR1dt = a1_prime + a2_prime - a3_prime + isqrV*(sqrt(a1_prime)*noise_a1(t) + sqrt(a2_prime)*noise_a2(t) - sqrt(a3_prime)*noise_a3(t));
        dR2dt = b1_prime + b2_prime - b3_prime + isqrV*(sqrt(b1_prime)*noise_b1(t) + sqrt(b2_prime)*noise_b2(t) - sqrt(b3_prime)*noise_b3(t));
        dR3dt = c1_prime + c2_prime - c3_prime + isqrV*(sqrt(c1_prime)*noise_c1(t) + sqrt(c2_prime)*noise_c2(t) - sqrt(c3_prime)*noise_c3(t));
        dR4dt = d1_prime + d2_prime - d3_prime + isqrV*(sqrt(d1_prime)*noise_d1(t) + sqrt(d2_prime)*noise_d2(t) - sqrt(d3_prime)*noise_d3(t));

        % apply state changes
        R1 += dR1dt * dt;
        R2 += dR2dt * dt;
        R3 += dR3dt * dt;
        R4 += dR4dt * dt;

    endfor

    figure; % plot each graphic on a separate window
    hold on;

    % scale data for easier visualization
    timescale = 1e5;
    quantscale = 1e-9;
    x = simulation / timescale;
    yR1 = R1s / quantscale;
    yR2 = R2s / quantscale;
    yR3 = R3s / quantscale;
    yR4 = R4s / quantscale;

    % plot results
    plot(x,yR1,'-m;R1;', x,yR2,'-k;R2;', x,yR3,'-r;R3;', x,yR4,'-g;R4;');
    xlabel("Time (10^5 seconds)");
    ylabel("Concentration (nM)");
    title("Stability Analysis (Reduced Model)");

    hold off;

    a = input("\nPress enter to exit ");
