tic;
%% Initialization
% close all; 
% clearvars;
% Setting paths
addpath(genpath(pwd));


%% Loading parameters
% fire_type = 'Syn';
% Ex_type = 'constant';
% params = model_parameters(fire_type);

for i = 1
    fire_type = 'Syn';
    Ex_type = 'constant';
    model_type = 'full_model';
    params = model_parameters(fire_type);
    params.Ex_Poisson_lambda = 7;
    params.M = 100;
    params.ne = 300;
    params.ni = 100;
    params.dt = 0.05;
    params.duration_time = 1000;
    params.tau_ee = 1.4;
    params.tau_ie = 1.2;
    params.tau_ei = 4.5;
    params.tau_i = 4.5;
    params.tau_r = 0;
    params.p_ee = 0.8; % P(preE activate postE)
    params.p_ie = 0.8; % P(preE activate postI)
    params.p_ei  = 0.8; % P(preI activate postE)
    params.p_ii = 0.8; % P(preI activate postI)
    %% 稳定点参数
    params.tau_ee = 4;
    params.tau_ie = 4;
    params.s_ee     = 1;
    params.s_ie     = 0.95;
    params.s_ei     = 2.61;
    params.s_ii     = 2.45;
    
    %% 吸引子参数
%     params.tau_ee = 1.4;
%     params.tau_ie = 1.2;
%     params.s_ee     = 0.94;
%     params.s_ie     = 1.25;
%     params.s_ei     = 2.61;
%     params.s_ii     = 2.45;


%     res_Ifull_model = run_ODEfull_model(params);
%     res_Ifull_model.Ex_type = Ex_type;
%     start.nfe = start.nfe+1;
    res_ODEfull_model = run_MIFIODE2_model(params);
    res_ODEfull_model.Ex_type = Ex_type;

%     save(['D:\matlab project\changexinhao\model_res\', strcat(fire_type,'_',Ex_type, model_type, 'See',,'var10exvar2_res.mat')], 'res_ODEfull_model', 'params')
%     clear('res_ODEfull_model');
%     clear('params');
end


toc