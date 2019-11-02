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
    simulation = 0 : dt : 1.5e6; % simulation length control

    % reaction initial conditions for each run
    SI =   1e-9 * [0.1, 0.1, 5 , 7.5, 10, 10]; % constant input
    SR12 = 1e-9 * [50 , 0  , 50, 50 , 50, 0 ];
    SR34 = 1e-9 * [0  , 50 , 0 , 0  , 0 , 50];


% ================================ MULTISIM ====================================

    % precomputed invariants
    alpha = Kt * Kb / Kd_m;
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


    % for each set of parameters, make a simulation run and print results
    for s = 1 : length(SI)

        % preallocate arrays to store concentrations over time
        Is = zeros(1, length(simulation));
        R1s = zeros(1, length(simulation));
        R2s = zeros(1, length(simulation));
        R3s = zeros(1, length(simulation));
        R4s = zeros(1, length(simulation));

        % initial state concentrations (M)
        Is(:) = SI(s); % constant input vector
        R1 = R2 = SR12(s);
        R3 = R4 = SR34(s);


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
            dR1dt = alpha * (HaI*HrR2 + HrR3) - Kd_p*R1;
            dR2dt = alpha * HrR4 * (HaI + HrR3) - Kd_p*R2;
            dR3dt = alpha * (HaI*HrR4 + HrR1) - Kd_p*R3;
            dR4dt = alpha * HrR2 * (HaI + HrR1) - Kd_p*R4;

            % apply state changes
            R1 += dR1dt * dt;
            R2 += dR2dt * dt;
            R3 += dR3dt * dt;
            R4 += dR4dt * dt;

        endfor


        figure; % plot each graphic in its separate window
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
        plot(simulation,R1s,'-m;R1;', simulation,R2s,'-k;R2;', simulation,R3s,'-r;R3;', simulation,R4s,'-g;R4;');
        xlabel("Time (10^5 seconds)");
        ylabel("Concentration (nM)");
        title("Stability Analysis");

        hold off;

    endfor
    a = input("\nPress enter to exit ");
