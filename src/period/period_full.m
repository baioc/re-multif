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
    simulation = 0 : dt : 8.0e5; % simulation length control

    % initial state concentrations (M)
    R1 = 50e-9;
    R2 = 50e-9;
    R3 = 0e-9;
    R4 = 0e-9;
    mR1P1 = 0e-9;
    mR4P1 = 0e-9;
    mR2P2 = 0e-9;
    mR3P2 = 0e-9;
    mR1P3 = 0e-9;
    mR2P4 = 0e-9;
    mR3P5 = 0e-9;
    mR4P6 = 0e-9;


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

    period_vector = zeros(1, 14);
    sin_vector = zeros(1, 14);

    var = -1;

    for s = 1:14

        low_value = 40e-9;
        high_value = 60e-9;
        is_period = 0;
        is_period_2 = 0;
        end_time_2 = -1;
        begin_time_2 = -1;
        time_stamp_begin = -1;
        time_stamp_end = -1;
        % input vector
        period = 0.18 + (s - 1)*0.12;
        period_time_2 = -1;
        if period == 0
            period_vector(s) = 0;
            continue
        endif
        Is = (6 .- 25*cos(simulation .* 2*pi/(period*1e5)) .+ 25) * 1e-9; % sinusoidal

        % actual simulation
        for t = 1 : length(simulation)

            % store current state
            R1s(t) = R1;
            R2s(t) = R2;
            R3s(t) = R3;
            R4s(t) = R4;

            if Is(t) > 50e-9
                if is_period_2 == 0
                    is_period_2 = 1;
                    if begin_time_2 == -1
                        begin_time_2 = t;
                    else
                        end_time_2 = t;
                        period_time_2 = (end_time_2 - begin_time_2) * 60;
                        begin_time_2 = -1;
                    endif
                endif
            endif
            if Is(t) < 15e-9 && is_period_2 == 1
                is_period_2 = 0;
            endif

            if R2 < low_value && is_period == 1
                is_period = 0;
                if time_stamp_begin == -1
                    time_stamp_begin = t;
                else
                    time_stamp_end = t;
                    h = (time_stamp_end - time_stamp_begin) * 60;
                    if h > period_vector(s)
                        period_vector(s) = h;
                    endif
                    time_stamp_begin = -1;
                endif
            endif

            if R2 > high_value && is_period == 0
                is_period = 1;
            endif

            % common subexpression optimization (used below)
            IsP1 = Is(t)^Na_P1;
            activationP1 = IsP1 / (Ka_P1 + IsP1);
            repressionP1 = 1 / (1 + (R2^Nr_R2P1 / Kr_R2P1));
            IsP2 = Is(t)^Na_P2;
            activationP2 = IsP2 / (Ka_P2 + IsP2);
            repressionP2 = 1 / (1 + (R4^Nr_R4P2 / Kr_R4P2));

            % compute variation
            dR1dt = Kt_R1*(mR1P1 + mR1P3) - Kd_R1*R1;
            dR2dt = Kt_R2*(mR2P2 + mR2P4) - Kd_R2*R2;
            dR3dt = Kt_R3*(mR3P2 + mR3P5) - Kd_R3*R3;
            dR4dt = Kt_R4*(mR4P1 + mR4P6) - Kd_R4*R4;
            dmR1P1dt = Kc_P1 + Kb_P1*activationP1*repressionP1 - Kd_mR1P1*mR1P1;
            dmR4P1dt = Kc_P1 + Kb_P1*activationP1*repressionP1 - Kd_mR4P1*mR4P1;
            dmR2P2dt = Kc_P2 + Kb_P2*activationP2*repressionP2 - Kd_mR2P2*mR2P2;
            dmR3P2dt = Kc_P2 + Kb_P2*activationP2*repressionP2 - Kd_mR3P2*mR3P2;
            dmR1P3dt = Kc_P3 + Kb_P3 * (1 / (1 + (R3^Nr_R3P3 / Kr_R3P3))) - Kd_mR1P3*mR1P3;
            dmR2P4dt = Kc_P4 + Kb_P4 * (1 / (1 + (R3^Nr_R3P4 / Kr_R3P4))) * (1 / (1 + (R4^Nr_R4P4 / Kr_R4P4))) - Kd_mR2P4*mR2P4;
            dmR3P5dt = Kc_P5 + Kb_P5 * (1 / (1 + (R1^Nr_R1P5 / Kr_R1P5))) - Kd_mR3P5*mR3P5;
            dmR4P6dt = Kc_P6 + Kb_P6 * (1 / (1 + (R1^Nr_R1P6 / Kr_R1P6))) * (1 / (1 + (R2^Nr_R2P6 / Kr_R2P6))) - Kd_mR4P6*mR4P6;

            % apply state changes
            R1 += dR1dt * dt;
            R2 += dR2dt * dt;
            R3 += dR3dt * dt;
            R4 += dR4dt * dt;
            mR1P1 += dmR1P1dt * dt;
            mR4P1 += dmR4P1dt * dt;
            mR2P2 += dmR2P2dt * dt;
            mR3P2 += dmR3P2dt * dt;
            mR1P3 += dmR1P3dt * dt;
            mR2P4 += dmR2P4dt * dt;
            mR3P5 += dmR3P5dt * dt;
            mR4P6 += dmR4P6dt * dt;

        endfor
        if period_vector(s) == 0
            period_vector(s) = period_time_2;
        endif
        sin_vector(s) = period_time_2;
    endfor


% ================================ RESULTS =====================================

    figure;
    hold on;

    %% scale data for easier visualization
    %timescale = 1e5;
    %quantscale = 1e-9;
    %x = simulation / timescale;
    %yI = Is / quantscale;
    %yR1 = R1s / quantscale;
    %yR2 = R2s / quantscale;
    %yR3 = R3s / quantscale;
    %yR4 = R4s / quantscale;

    % model subplot
    plot(sin_vector, period_vector, '@r');
    xlabel("Sin(x) period increment");
    ylabel("Period increase");
    title("Full Model");

    hold off;
    print('period.png');
    a = input("\nPress enter to exit ");