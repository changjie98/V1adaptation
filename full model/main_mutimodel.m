tic;
%% Initialization
% close all; 
clearvars;
% Setting paths
addpath(genpath(pwd));


%% Loading parameters
% fire_type = 'Syn';
% Ex_type = 'constant';
% params = model_parameters(fire_type);


for i = 2
    fire_type = 'Syn';
    Ex_type = 'constant';
    model_type = 'full_model';
    params = model_parameters(fire_type);
    params.Ex_Poisson_lambda = 7;
    params.ne = 75;
    params.ni = 25;
    params.duration_time = 1000;
%     params.tau_ee = 4;
    params.tau_r = 2;
%     params.p_ee = 0.15; % P(preE activate postE)
%     params.p_ie = 0.5; % P(preE activate postI)
%     params.p_ei  = 0.5; % P(preI activate postE)
%     params.p_ii = 0.4; % P(preI activate postI)
%     params.s_ee     = 5;
%     params.s_ie     = 2;
%     params.s_ei     = 4.91;
%     params.s_ii     = 4.91;
    res_full_model = run_adpfull_model(params, 'constant');
%     res_full_model.firing_rate = firing_rate(res_full_model,params,'part');
    res_full_model.Ex_type = Ex_type;
%     save(['D:\matlab project\chengexinhao\model_res\limitcircle\N\200\', strcat(fire_type,'_',Ex_type, model_type, 'networksize200See',num2str(i),'_res.mat')], 'res_full_model', 'params')
%     clear('res_full_model');
%     clear('params');
end


% for i = 1:4
%     fire_type = 'Syn';
%     Ex_type = 'constant';
%     model_type = 'full_model';
%     params = model_parameters(fire_type);
%     params.Ex_Poisson_lambda = 7;
%     params.ne = 300;
%     params.ni = 100;
%     params.duration_time = 5000;
%     
%     params.tau_r = 2;
%     params.p_ee = 0.15; 
%     params.p_ie = 0.5; 
%     params.p_ei  = 0.5; 
%     params.p_ii = 0.4;
%     params.s_ee     = 5;
%     params.s_ie     = i;
%     params.s_ei     = 4.91;
%     params.s_ii     = 4.91;
%     res_full_model = run_full_model(params, 'constant');
%     firing_rate = firing_rate(res_full_model, params);
%     res_full_model.firing_rate = firing_rate;
%     res_full_model.Ex_type = Ex_type;
%     save(['D:\matlab project\chengexinhao\model_res\limitcircle\sie\', strcat(fire_type,'_',Ex_type, model_type, 'networksize400_taur2sie',num2str(i),'_reslongtime.mat')], 'res_full_model', 'params')
%     clear('res_full_model');
%     clear('params');
%     clear('firing_rate');
% end
%% run constant full models
% for i = 1:4
%     fire_type = 'Syn';
%     Ex_type = 'constant';
%     model_type = 'full_model';
%     params = model_parameters(fire_type);
%     params.Ex_Poisson_lambda = i;
%     params.ne = 300;
%     params.ni = 100;
%     params.duration_time = 5000;
% %     params.tau_r = 5;
%     params.p_ee = 0.15; 
%     params.p_ie = 0.5; 
%     params.p_ei  = 0.2+0.1*i; 
%     params.p_ii = 0.4;
%     params.s_ee     = 5;
%     params.s_ie     = 2;
%     params.s_ei     = 4.91;
%     params.s_ii     = 4.91;
%     res_full_model = run_full_model(params, 'constant');
%     firing_rate = firing_rate(res_full_model, params);
%     res_full_model.firing_rate = firing_rate;
%     res_full_model.Ex_type = Ex_type;
%     save(['D:\matlab project\chengexinhao\model_res\limitcircle\pei\', strcat(fire_type,'_',Ex_type, model_type, 'networksize400_pei',num2str(0.2+0.1*i),'_reslongtime.mat')], 'res_full_model', 'params')
%     clear('res_full_model');
%     clear('params');
%     clear('firing_rate');
% end

% for i = 7
%     fire_type = 'Syn';
%     Ex_type = 'constant';
%     model_type = 'full_model';
%     params = model_parameters(fire_type);
%     params.Ex_Poisson_lambda = i;
%     params.tau_ee = 1.2;
%     params.tau_ie = 1.4;
%     params.tau_i = 4.5;
%     params.ne = 300;
%     params.ni = 100;
%     params.duration_time = 5000;
% %     params.tau_r = 5;
% %     params.p_ee = 0.2; 
% %     params.p_ie = 0.5; 
% %     params.p_ei  = 0.5; 
% %     params.p_ii = 0.4;
% %     params.M = 50*i;
%     params.s_ee     = 4;
%     params.s_ie     = 2;
%     params.s_ei     = 4.91;
%     params.s_ii     = 4.91;
%     res_full_model = run_full_model(params, 'constant');
%     firing_rate = firing_rate(res_full_model, params);
%     res_full_model.firing_rate = firing_rate;
%     res_full_model.Ex_type = Ex_type;
%     save(['D:\matlab project\chengexinhao\model_res\limitcircle\', strcat(fire_type,'_',Ex_type, model_type, 'networksize400_see5ex',num2str(i),'_reslongtime.mat')], 'res_full_model', 'params')
%     clear('res_full_model');
%     clear('params');
%     clear('firing_rate');
% end

%% run constant reduced models
% for i = 2:12
%     fire_type = 'Syn';
%     Ex_type = 'constant';
%     model_type = 'reduced_model';
%     params = model_parameters(fire_type);
%     params.Ex_Poisson_lambda = i;
%     res_reduced_model = run_reduced_model(params, 'constant');
%     firing_rate = firing_rate(res_reduced_model, params);
%     res_reduced_model.firing_rate = firing_rate;
%     save(['D:\matlab\bin\change xinhao\model_res\', strcat(fire_type,'_',Ex_type,num2str(i),'model_type','_res.mat')], 'res_reduced_model', 'params')
%     clearvars;
% end


%% run constant CG models
% for i = 2:12
%     fire_type = 'Syn';
%     Ex_type = 'constant';
%     model_type = 'reduced_model';
%     params = model_parameters(fire_type);
%     params.S_ee = 11.25;
%     params.S_ie = 12.5;
%     params.S_e = 24;
%     params.S_i = 48;
%     params.S_ii = 10;
%     params.S_ei = 38;
%     params.a_ee=0.5;
%     params.a_ie=0.5;
%     params.a_ei=0.79;
%     params.a_ii=0.21;
%     params.Ex_Poisson_lambda = i;
%     res_CG_model = run_CG_model(params, 'constant');
%     firing_rate = firing_rate(res_CG_model, params);
%     res_CG_model.firing_rate = firing_rate;
%     save(['D:\matlab\bin\change xinhao\model_res\', strcat(fire_type,'_',Ex_type,num2str(i),'model_type','_res.mat')], 'res_CG_model', 'params')
%     clearvars;
% end

%% run sin full models

% for i = 7:7
%     fire_type = 'Syn';
%     Ex_type = 'sin';
%     model_type = 'full_model';
%     params = model_parameters(fire_type);
%     C2 = [0.001*pi, 0.01*pi, 0.03*pi, 0.05*pi, 0.07*pi, 0.09*pi, 0.1*pi, 0.11*pi, 0.13*pi, 0.15*pi, 0.2*pi, 0.5*pi, pi];
%     params.Ex_sin_C1 = 5; 
%     params.Ex_sin_C2 = C2(i);
%     params.Ex_sin_C3 = 0;
%     params.Ex_sin_C4 = 7;
%     params.duration_time = 5000;
%     
%     params.ne = 300;
%     params.ni = 100;
%     params.tau_r = 2;
% %     params.tau_ee = 4;
% %     params.p_ee = 0.075; 
% %     params.p_ie = 0.25; 
% %     params.p_ei  = 0.25; 
% %     params.p_ii = 0.2;
%     params.s_ee     = 6;
%     params.s_ie     = 2;
%     params.s_ei     = 4.91;
%     params.s_ii     = 4.91;
%     res_full_model = run_full_model(params, Ex_type);
%     firing_rate = firing_rate(res_full_model, params);
%     res_full_model.firing_rate = firing_rate;
%     res_full_model.Ex_type = 'sin';
%     save(['D:\matlab project\chengexinhao\model_res\limitcircle\adp\', strcat(fire_type,'_',Ex_type, model_type, 'sinSyn',num2str(400),'_res.mat')], 'res_full_model', 'params')
%     clear('res_full_model');
%     clear('params');
%     clear('firing_rate');
% end



%% run sin reduced models
% for i = 2:12
%     fire_type = 'Syn';
%     Ex_type = 'sin';
%     model_type = 'reduced_model';
%     load('nonref_P.mat')
%     params = model_parameters(fire_type);
%     C2 = [0.001*pi, 0.01*pi, 0.03*pi, 0.05*pi, 0.07*pi, 0.09*pi, 0.1*pi, 0.11*pi, 0.13*pi, 0.15*pi, 0.2*pi, 0.5*pi, pi];
%     params.Ex_sin_C2=C2(i);
%     res_reduced_model = run_reduced_model(params, Ex_type, P3_stat_noref);
%     firing_rate = firing_rate(res_reduced_model, params);
%     res_reduced_model.firing_rate = firing_rate;
%     save(['D:\matlab\bin\change xinhao\model_res\', strcat(fire_type,'_',Ex_type,num2str(i),'model_type','_res.mat')], 'res_reduced_model', 'params')
%     clearvars;
% end


%% run sin CG models
% for i = 1:13
%     fire_type = 'Syn';
%     Ex_type = 'sin';
%     model_type = 'CG_model';
%     params = model_parameters(fire_type);
%     params.S_ee = 11.25;
%     params.S_ie = 12.5;
%     params.S_e = 24;
%     params.S_i = 48;
%     params.S_ii = 10;
%     params.S_ei = 38;
%     params.a_ee=0.5;
%     params.a_ie=0.5;
%     params.a_ei=0.79;
%     params.a_ii=0.21;
%     params.Ex_Poisson_lambda = i;
%     C2 = [0.001*pi, 0.01*pi, 0.03*pi, 0.05*pi, 0.07*pi, 0.09*pi, 0.1*pi, 0.11*pi, 0.13*pi, 0.15*pi, 0.2*pi, 0.5*pi, pi];
%     params.Ex_sin_C2=C2(i);
%     load('nonref_P.mat');
%     res_CG_model = run_CG_model(params, 'sin', P3_stat_noref);
%     firing_rate = firing_rate(res_CG_model, params);
%     res_CG_model.firing_rate = firing_rate;
%     save(['D:\matlab project\chengexinhao\model_res\', strcat(fire_type,'_',Ex_type,num2str(C2(i)), '_', model_type, '_res.mat')], 'res_CG_model', 'params')
%     clearvars;
% end


toc