#!/usr/bin/octave -qf
% usage: octave save_parameters.m # then rename generated "parameters.mat" file

% experimental context parameters (Stochastic only works w/ single || switch)
experiment_class = 'single'; % 'single' | 'period' | 'oscillator' | 'switch'
model = 'Stochastic'; % 'QSSA' | 'Full' | 'Stochastic'
method = 'Euler'; % 'Euler' | 'ode45'

% only used in the 'Stochastic' model
random_seed =  0.8977812586558097;
noise_scaling = 1/26;

% only used in 'period' experiments
period_experiment_range = [2e4, 0.5e4, 16e4]; % min : step : max

% only used in 'oscillator' experiments
oscillator_input_range = [0.4e-9, 0.1e-9, 7e-9]; % min : step : max

% only used in 'switch' experiments
switch_trigger_delta = +39994e-10; % k_r
switch_trigger_time = [1.50e5 1.55e5]; % (s) start, stop
switch_R = [1 1 0 0]; % affects R1? R2? R3? R4?

% initial concentrations of repressor proteins
R1 = 50e-9; % (M)
R2 = 50e-9; % (M)
R3 = 0e-9; % (M)
R4 = 0e-9; % (M)

% input signal configuration
input_type = 'sine'; % 'sine' | 'square' | 'constant'

% when input is "constant", this is its value (except in 'oscillator')
input_dc_level = 6e-9; % (M)

% these are only used in 'sine' and 'square' inputs
input_period = 0.9e5; % (s), will be overrided in 'period' experiments
input_amplitude = 50e-9; % (M)
input_duty_cycle = 0.5; % (normalized %)

% simulation time controls
simulation_total_time = 4.0e5; %(s), overrided in all but 'single' | 'switch'
simulation_step = 60; %(s)

% plot configuration
plot_output_filename = "stochastic-freqdiv.pdf"; % saved plot filename
plot_include_input_signal = false; % (bool), only used in 'single' experiments
plot_aspect = [1 0.334 1]; % vector or "auto" that is passed to pbaspect()

% save experiment configuration in plain text
save "parameters.mat" ...
     experiment_class model method ...
     random_seed noise_scaling ...
     period_experiment_range ...
     oscillator_input_range ...
     switch_trigger_delta switch_trigger_time switch_R ...
     R1 R2 R3 R4 ...
     input_type input_dc_level input_period input_amplitude input_duty_cycle ...
     simulation_total_time simulation_step ...
     plot_output_filename plot_include_input_signal plot_aspect;
