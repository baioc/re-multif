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

% activation constants (M)
Ka_P1 = 2e-8;
Ka_P2 = 2e-8;

% repression constants (M)
Kr_R2P1 = 6e-10;
Kr_R4P2 = 6e-10;
Kr_R3P3 = 6e-10;
Kr_R3P4 = 6e-10;
Kr_R4P4 = 6e-10;
Kr_R1P5 = 6e-10;
Kr_R1P6 = 6e-10;
Kr_R2P6 = 6e-10;

% maximum/unrepressed transcription rates (M/s)
Kb_P1 = 4e-10;
Kb_P2 = 4e-10;
Kb_P3 = 4e-10;
Kb_P4 = 4e-10;
Kb_P5 = 4e-10;
Kb_P6 = 4e-10;

% simplified coefficients
a1 = Kt_R1 * Kb_P1 / Kd_mR1P1;
a2 = Kt_R1 * Kb_P3 / Kd_mR1P3;
b1 = Kt_R2 * Kb_P2 / Kd_mR2P2;
b2 = Kt_R2 * Kb_P4 / Kd_mR2P4;
c1 = Kt_R3 * Kb_P2 / Kd_mR3P2;
c2 = Kt_R3 * Kb_P5 / Kd_mR3P5;
d1 = Kt_R4 * Kb_P1 / Kd_mR4P1;
d2 = Kt_R4 * Kb_P6 / Kd_mR4P6;

% abstractions for activation & repression Hill-function
Ha = @(X, n, K) X^n / (K + X^n);
Hr = @(X, n, K) 1 / (1 + (X^n / K));

% simulation parameters
dt = 1;
simulation = 0 : dt : 3.5e5;

% system initial state (M)
R1 = 50e-9;
R2 = 50e-9;
R3 = 0e-9;
R4 = 0e-9;

% concentrations over time
Is = (6 .- 25*cos(simulation .* 2*pi/0.9e5) .+ 25) * 1e-9;
R1s = zeros(1, length(simulation));
R2s = zeros(1, length(simulation));
R3s = zeros(1, length(simulation));
R4s = zeros(1, length(simulation));

% simulation
for t = 1 : length(simulation)

    R1s(t) = R1;
    R2s(t) = R2;
    R3s(t) = R3;
    R4s(t) = R4;

    dR1dt = a1*Ha(Is(t), Na_P1, Ka_P1)*Hr(R2, Nr_R2P1, Kr_R2P1) + a2*Hr(R3, Nr_R3P3, Kr_R3P3) - Kd_R1*R1;
    dR2dt = b1*Ha(Is(t), Na_P2, Ka_P2)*Hr(R4, Nr_R4P2, Kr_R4P2) + b2*Hr(R3, Nr_R3P4, Kr_R3P4)*Hr(R4, Nr_R4P4, Kr_R4P4) - Kd_R2*R2;
    dR3dt = c1*Ha(Is(t), Na_P2, Ka_P2)*Hr(R4, Nr_R4P2, Kr_R4P2) + c2*Hr(R1, Nr_R1P5, Kr_R1P5) - Kd_R3*R3;
    dR4dt = d1*Ha(Is(t), Na_P1, Ka_P1)*Hr(R2, Nr_R2P1, Kr_R2P1) + d2*Hr(R1, Nr_R1P6, Kr_R1P6)*Hr(R2, Nr_R2P6, Kr_R2P6) - Kd_R4*R4;

    R1 += dR1dt * dt;
    R2 += dR2dt * dt;
    R3 += dR3dt * dt;
    R4 += dR4dt * dt;

endfor

% plotting results
figure;
hold on;

    title("Reduced Model");

    subplot(2,1, 1);
    plot(simulation,R1s,'-m;R1;', simulation,R2s,'-k;R2;', simulation,R3s,'-r;R3;', simulation,R4s,'-g;R4;');
    xlabel("Time (s)");
    ylabel("Concentration (M)");

    subplot(2,1, 2);
    plot(simulation, Is, 'b;I;');
    xlabel("Time (s)");
    ylabel("Concentration (M)");

hold off;
a = input("\nPress enter to exit ");
