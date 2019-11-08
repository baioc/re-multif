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
    simulation = 0 : dt : 1.4e6; % simulation length control

    % input vector
    Is = zeros(1, length(simulation));

    % reaction initial conditions for each run
    names = {"1a", "1b", "2", "3", "4c_1", "4c_2"};
    multiI   = 1e-9 * [0.1, 0.1, 5 , 7.5, 10, 10]; % constant input concentration
    multiR12 = 1e-9 * [50 , 0  , 50, 50 , 50, 0 ];
    multiR34 = 1e-9 * [0  , 50 , 0 , 0  , 0 , 50];


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
    for sim = 1 : length(multiI)

        % set initial state concentrations
        Is(:) = multiI(sim);
        R1 = multiR12(sim);
        R2 = multiR12(sim);
        R3 = multiR34(sim);
        R4 = multiR34(sim);

        % preallocate arrays to store concentrations over time
        R1s = zeros(1, length(simulation));
        R2s = zeros(1, length(simulation));
        R3s = zeros(1, length(simulation));
        R4s = zeros(1, length(simulation));

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

            % compute variation
            dR1dt = Kc_P1+Kc_P3 + a1*activationP1*repressionP1 + a2*Hr(R3, Nr_R3P3, Kr_R3P3) - Kd_R1*R1;
            dR2dt = Kc_P2+Kc_P4 + b1*activationP2*repressionP2 + b2*Hr(R3, Nr_R3P4, Kr_R3P4)*Hr(R4, Nr_R4P4, Kr_R4P4) - Kd_R2*R2;
            dR3dt = Kc_P2+Kc_P5 + c1*activationP2*repressionP2 + c2*Hr(R1, Nr_R1P5, Kr_R1P5) - Kd_R3*R3;
            dR4dt = Kc_P1+Kc_P6 + d1*activationP1*repressionP1 + d2*Hr(R1, Nr_R1P6, Kr_R1P6)*Hr(R2, Nr_R2P6, Kr_R2P6) - Kd_R4*R4;

            % apply state changes
            R1 += dR1dt * dt;
            R2 += dR2dt * dt;
            R3 += dR3dt * dt;
            R4 += dR4dt * dt;

        endfor

        figure; % plot each graphic separately
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
        title(cstrcat("Bifurcation Analysis experiment ", names{sim})); %, " (Reduced Model)"));

        hold off;
        names{sim}
        maxR1 = max(yR1(6e5/dt:end))
        maxR4 = max(yR4(6e5/dt:end))
        print(strcat('bifurcation-', names{sim}, '-qssa.pdf')); % put in the folder the script is run from

    endfor
