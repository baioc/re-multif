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
    dt = 60;                   % 1 minute timesteps
    simulation = 0 : dt : 4e5; % simulation length control

    % constant input vector
    Is = zeros(1, length(simulation));
    Is(:) = 50e-9;

    % induced decrease in repressor binding affinity used to trigger the switch
    switch_trigger = 4e-6 - 6e-10;

    % save original Hill repression constants
    Kr_R2P1_ori = Kr_R2P1;
    Kr_R4P2_ori = Kr_R4P2;
    Kr_R3P3_ori = Kr_R3P3;
    Kr_R3P4_ori = Kr_R3P4;
    Kr_R4P4_ori = Kr_R4P4;
    Kr_R1P5_ori = Kr_R1P5;
    Kr_R1P6_ori = Kr_R1P6;
    Kr_R2P6_ori = Kr_R2P6;

    % initial state concentrations (M)
    R1 = R2 = 50e-9;
    R3 = R4 = 0e-9;
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
    T12s = zeros(1, length(simulation));
    T34s = zeros(1, length(simulation));

    % actual simulation
    for t = 1 : length(simulation)

        % store current state
        R1s(t) = R1;
        R2s(t) = R2;
        R3s(t) = R3;
        R4s(t) = R4;

        % change binding following hold or release of switch trigger
        if t >= 1e5 / dt && t <= 1.05e5 / dt
            Kr_R1P5 = Kr_R1P5_ori + switch_trigger;
            Kr_R1P6 = Kr_R1P6_ori + switch_trigger;
            Kr_R2P6 = Kr_R2P6_ori + switch_trigger;
            Kr_R2P1 = Kr_R2P1_ori + switch_trigger;
            T12s(t) = Kr_R2P1_ori + switch_trigger;
            Kr_R3P3 = Kr_R3P3_ori;
            Kr_R3P4 = Kr_R3P4_ori;
            Kr_R4P4 = Kr_R4P4_ori;
            Kr_R4P2 = Kr_R4P2_ori;
        elseif t >= 2.50e5 / dt && t <= 2.55e5 / dt
            Kr_R1P5 = Kr_R1P5_ori;
            Kr_R1P6 = Kr_R1P6_ori;
            Kr_R2P6 = Kr_R2P6_ori;
            Kr_R2P1 = Kr_R2P1_ori;
            Kr_R3P3 = Kr_R3P3_ori + switch_trigger;
            Kr_R3P4 = Kr_R3P4_ori + switch_trigger;
            Kr_R4P4 = Kr_R4P4_ori + switch_trigger;
            Kr_R4P2 = Kr_R4P2_ori + switch_trigger;
            T34s(t) = Kr_R4P2_ori + switch_trigger;
        else
            Kr_R2P1 = Kr_R2P1_ori;
            Kr_R4P2 = Kr_R4P2_ori;
            Kr_R3P3 = Kr_R3P3_ori;
            Kr_R3P4 = Kr_R3P4_ori;
            Kr_R4P4 = Kr_R4P4_ori;
            Kr_R1P5 = Kr_R1P5_ori;
            Kr_R1P6 = Kr_R1P6_ori;
            Kr_R2P6 = Kr_R2P6_ori;
            T12s(t) = Kr_R2P1_ori;
            T34s(t) = Kr_R4P2_ori;
        endif

        % common subexpression optimization (used below)
        IsP1 = Is(t)^Na_P1;
        activationP1 = IsP1 / (Ka_P1^Na_P1 + IsP1);
        repressionP1 = 1 / (1 + (R2/Kr_R2P1)^Nr_R2P1);
        IsP2 = Is(t)^Na_P2;
        activationP2 = IsP2 / (Ka_P2^Na_P2 + IsP2);
        repressionP2 = 1 / (1 + (R4/Kr_R4P2)^Nr_R4P2);

        % compute variation
        dR1dt = Kt_R1*(mR1P1 + mR1P3) - Kd_R1*R1;
        dR2dt = Kt_R2*(mR2P2 + mR2P4) - Kd_R2*R2;
        dR3dt = Kt_R3*(mR3P2 + mR3P5) - Kd_R3*R3;
        dR4dt = Kt_R4*(mR4P1 + mR4P6) - Kd_R4*R4;
        dmR1P1dt = Kc_P1 + Kb_P1*activationP1*repressionP1 - Kd_mR1P1*mR1P1;
        dmR4P1dt = Kc_P1 + Kb_P1*activationP1*repressionP1 - Kd_mR4P1*mR4P1;
        dmR2P2dt = Kc_P2 + Kb_P2*activationP2*repressionP2 - Kd_mR2P2*mR2P2;
        dmR3P2dt = Kc_P2 + Kb_P2*activationP2*repressionP2 - Kd_mR3P2*mR3P2;
        dmR1P3dt = Kc_P3 + Kb_P3 * (1 / (1 + (R3/Kr_R3P3)^Nr_R3P3)) - Kd_mR1P3*mR1P3;
        dmR2P4dt = Kc_P4 + Kb_P4 * (1 / (1 + (R3/Kr_R3P4)^Nr_R3P4)) * (1 / (1 + (R4/Kr_R4P4)^Nr_R4P4)) - Kd_mR2P4*mR2P4;
        dmR3P5dt = Kc_P5 + Kb_P5 * (1 / (1 + (R1/Kr_R1P5)^Nr_R1P5)) - Kd_mR3P5*mR3P5;
        dmR4P6dt = Kc_P6 + Kb_P6 * (1 / (1 + (R1/Kr_R1P6)^Nr_R1P6)) * (1 / (1 + (R2/Kr_R2P6)^Nr_R2P6)) - Kd_mR4P6*mR4P6;

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


% ================================ RESULTS =====================================

    figure;
    hold on;

    % scale data for easier visualization
    timescale = 1e5;
    quantscale = 1e-9;
    x = simulation / timescale;
    yR1 = R1s / quantscale;
    yR2 = R2s / quantscale;
    yR3 = R3s / quantscale;
    yR4 = R4s / quantscale;

    % result subplot
    subplot(2, 1, 1);
    plot(x,yR1,'--m;R1;', x,yR2,':k;R2;', x,yR3,'-r;R3;', x,yR4,'-.g;R4;');
    pbaspect([1 0.334 1]);
    legend('location', 'east');
    ylabel("Concentration (nM)");

    % trigger subplot
    subplot(2, 1, 2);
    plot(x,T12s,'-b;k_A for R1 & R2;', x,T34s,'-r;k_A for R3 & R4;');
    axis([-Inf, +Inf, 0, max(T12s)*1.1]);
    xlabel("Time (10^5 seconds)");
    ylabel("Repressor k_A (M)");

    hold off;
    print(strcat('switch-full-', num2str(Is(1)/quantscale), '_.pdf')); % put in the folder the script is run from
