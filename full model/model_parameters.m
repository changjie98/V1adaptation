function [params] = model_parameters(need_set)
params = struct;


%% Neuron
params.ne = 75; % Number of E neurons
params.ni = 25; % Number of I neurons
params.M = 100; % threshold
params.Mr = 66; 
params.tau_r = 2;
params.tau_adp = 100;

%% Timeline
params.duration_time = 1000;
params.time_delta = 1;


%% Ex input
params.Ex_sin_C1 = 5; 
params.Ex_sin_C2 = 0.01;
params.Ex_sin_C3 = 0;
params.Ex_sin_C4 = 7;
params.Ex_Poisson_lambda = 7;


%% Neuron connections
params.p_ee = 0.15; % P(preE activate postE)
params.p_ie = 0.5; % P(preE activate postI)
params.p_ei  = 0.5; % P(preI activate postE)
params.p_ii = 0.4; % P(preI activate postI)
params.s_ee  = 25; % synaptic strength
params.s_ie = 10;
params.s_ei  = 25;
params.s_ii  = 25;

switch need_set 
    case{'Hom'}  
    params.tau_ee = 4;
    params.tau_ie = 1.2;
    params.tau_i = 4.5;
    case{'Reg'}
    params.tau_ee = 2;
    params.tau_ie = 1.2;
    params.tau_i = 4.5;
    case{'Syn'}
    params.tau_ee = 1.4;
    params.tau_ie = 1.2;
    params.tau_i = 4.5;
end

%% performance
% cluster; Spiking Volley Detection, looking for MFE
params.cluster_delta = 0.33;
params.cluster_eps   = 8;

% SSI: Spike Synchrony Index, Time to calculate the average intensity of the spike
params. w = 5;

% spectrogram: power spectrum density, How does the power of the signal vary with frequency
params.sdbin = 2.5;
params.spectrogram_timewindow = 100;
params.frequency_range = [20,80];


end

