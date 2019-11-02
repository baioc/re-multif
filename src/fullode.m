clear;
format long;

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

% Hill coefficients
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

% Hill activation constants (M^n)
Ka_P1 = 2e-8 ^ Na_P1;
Ka_P2 = 2e-8 ^ Na_P2;

% Hill repression constants (M^n)
Kr_R2P1 = 6e-10 ^ Nr_R2P1;
Kr_R4P2 = 6e-10 ^ Nr_R4P2;
Kr_R3P3 = 6e-10 ^ Nr_R3P3;
Kr_R3P4 = 6e-10 ^ Nr_R3P4;
Kr_R4P4 = 6e-10 ^ Nr_R4P4;
Kr_R1P5 = 6e-10 ^ Nr_R1P5;
Kr_R1P6 = 6e-10 ^ Nr_R1P6;
Kr_R2P6 = 6e-10 ^ Nr_R2P6;

% maximum/unrepressed transcription rates (M/s)
Kb_P1 = 4e-10;
Kb_P2 = 4e-10;
Kb_P3 = 4e-10;
Kb_P4 = 4e-10;
Kb_P5 = 4e-10;
Kb_P6 = 4e-10;

% simulation parameters
dt = 1;
simulation = 0 : dt : 4.0e5;

% system initial state (M)
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

% concentrations over time
Is = (6 .- 25*cos(simulation .* 2*pi/0.9e5) .+ 25) * 1e-9;
R1s = zeros(1, length(simulation));
R2s = zeros(1, length(simulation));
R3s = zeros(1, length(simulation));
R4s = zeros(1, length(simulation));
mR1P1s = zeros(1, length(simulation));
mR4P1s = zeros(1, length(simulation));
mR2P2s = zeros(1, length(simulation));
mR3P2s = zeros(1, length(simulation));
mR1P3s = zeros(1, length(simulation));
mR2P4s = zeros(1, length(simulation));
mR3P5s = zeros(1, length(simulation));
mR4P6s = zeros(1, length(simulation));

% simulation
for t = 1 : length(simulation)

    R1s(t) = R1;
    R2s(t) = R2;
    R3s(t) = R3;
    R4s(t) = R4;
    mR1P1s(t) = mR1P1;
    mR4P1s(t) = mR4P1;
    mR2P2s(t) = mR2P2;
    mR3P2s(t) = mR3P2;
    mR1P3s(t) = mR1P3;
    mR2P4s(t) = mR2P4;
    mR3P5s(t) = mR3P5;
    mR4P6s(t) = mR4P6;

    dR1dt = Kt_R1*(mR1P1 + mR1P3) - Kd_R1*R1;
    dR2dt = Kt_R2*(mR2P2 + mR2P4) - Kd_R2*R2;
    dR3dt = Kt_R3*(mR3P2 + mR3P5) - Kd_R3*R3;
    dR4dt = Kt_R4*(mR4P1 + mR4P6) - Kd_R4*R4;
    dmR1P1dt = Kb_P1 * (Is(t)^Na_P1 / (Ka_P1 + Is(t)^Na_P1)) * (1 / (1 + (R2^Nr_R2P1 / Kr_R2P1))) - Kd_mR1P1*mR1P1;
    dmR4P1dt = Kb_P1 * (Is(t)^Na_P1 / (Ka_P1 + Is(t)^Na_P1)) * (1 / (1 + (R2^Nr_R2P1 / Kr_R2P1))) - Kd_mR4P1*mR4P1;
    dmR2P2dt = Kb_P2 * (Is(t)^Na_P2 / (Ka_P2 + Is(t)^Na_P2)) * (1 / (1 + (R4^Nr_R4P2 / Kr_R4P2))) - Kd_mR2P2*mR2P2;
    dmR3P2dt = Kb_P2 * (Is(t)^Na_P2 / (Ka_P2 + Is(t)^Na_P2)) * (1 / (1 + (R4^Nr_R4P2 / Kr_R4P2))) - Kd_mR3P2*mR3P2;
    dmR1P3dt = Kb_P3 * (1 / (1 + (R3^Nr_R3P3 / Kr_R3P3))) - Kd_mR1P3*mR1P3;
    dmR2P4dt = Kb_P4 * (1 / (1 + (R3^Nr_R3P4 / Kr_R3P4))) * (1 / (1 + (R4^Nr_R4P4 / Kr_R4P4))) - Kd_mR2P4*mR2P4;
    dmR3P5dt = Kb_P5 * (1 / (1 + (R1^Nr_R1P5 / Kr_R1P5))) - Kd_mR3P5*mR3P5;
    dmR4P6dt = Kb_P6 * (1 / (1 + (R1^Nr_R1P6 / Kr_R1P6))) * (1 / (1 + (R2^Nr_R2P6 / Kr_R2P6))) - Kd_mR4P6*mR4P6;

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

% plotting results
figure;
hold on;

    timescale = 1e5;
    quantscale = 1e-9;

    simulation /= timescale;
    Is  /= quantscale;
    R1s /= quantscale;
    R2s /= quantscale;
    R3s /= quantscale;
    R4s /= quantscale;

    subplot(2,1, 1);
    plot(simulation,R1s,'-m;R1;', simulation,R2s,'-k;R2;', simulation,R3s,'-r;R3;', simulation,R4s,'-g;R4;');
    xlabel("Time (10^5 seconds)");
    ylabel("Concentration (nM)");
    title("Full Model");

    subplot(2,1, 2);
    plot(simulation, Is, 'b;I;');
    xlabel("Time (10^5 seconds)");
    ylabel("Concentration (nM)");

hold off;
a = input("\nPress enter to exit ");
