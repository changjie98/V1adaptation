tic
%% 模型包括
%% 设置参数
params.J_ex = 7;
params.M = 100;
params.Mr = 66;
params.ne = 300;
params.ni = 100;
params.dt = 0.01;
params.duration_time = 500;
params.tau_ee = 1.4;
params.tau_ie = 1.2;
params.tau_i = 4.5;
params.tau_r = 0;
params.p_ee = 0.4; % P(preE activate postE)
params.p_ie = 0.4; % P(preE activate postI)
params.p_ei = 0.4; % P(preI activate postE)
params.p_ii = 0.4; % P(preI activate postI)
params.s_ee = 0.95;
params.s_ie = 2.8;
params.s_ei = 2.2;
params.s_ii = 2;
params.V_bin = 10;
params.V_bin_min = -4;
params.V_bin_num = params.M /params.V_bin - params.V_bin_min;
params.digit_num = 10;
params.t_end = params.duration_time/params.dt;

%% 初始化几个变量
% 每个res包含如下变量：V all bin 、Neuron number in each bin、fr、ref、H（mean、var）、I（mean、var）
V_e_all = zeros(params.t_end, params.V_bin_num);
V_e_mean = zeros(params.t_end, params.V_bin_num);
V_e_all(1,5) = 0.0001; % 默认没有给初始状态的话，就设定为神经元电压都为0.0001
n_e = zeros(params.t_end, params.V_bin_num);
n_e(1,5) = params.ne;
V_e_mean(1,:) = V_e_all(1,:)./n_e(1,:);
fr_e = zeros(params.t_end, 1);
ref_e = zeros(params.t_end, 1);
H_ee_mean = zeros(params.t_end, 1);
H_ei_mean = zeros(params.t_end, 1);
H_ee_var = zeros(params.t_end, 1);
H_ei_var = zeros(params.t_end, 1);
I_e_mean = zeros(params.t_end, params.V_bin_num);
I_e_var = zeros(params.t_end, params.V_bin_num);

V_i_all = zeros(params.t_end, params.V_bin_num);
V_i_mean = zeros(params.t_end, params.V_bin_num);
V_i_all(1,5) = 0.0001;
n_i = zeros(params.t_end, params.V_bin_num);
n_i(1,5) = params.ni;
V_i_all = V_i_all./n_i;
V_i_mean(1,:) = V_i_all(1,:)./n_i(1,:);
fr_i = zeros(params.t_end, 1);
ref_i = zeros(params.t_end, 1);
H_ie_mean = zeros(params.t_end, 1);
H_ii_mean = zeros(params.t_end, 1);
H_ie_var = zeros(params.t_end, 1);
H_ii_var = zeros(params.t_end, 1);
I_i_mean = zeros(params.t_end, params.V_bin_num);
I_i_var = zeros(params.t_end, params.V_bin_num);

%% 可以调整输入为各种形式，这里用了最简单的稳定泊松输入，均值等于方差
I_eex = zeros(params.t_end,params.V_bin_num) + params.J_ex;
I_iex = zeros(params.t_end,params.V_bin_num) + params.J_ex;

%% 配套get_start.m使用，如果有特殊的初始状态，调整havestart=1
havestart = 1;
if havestart
    H_ee_mean(1) = start.H_ee_mean;
    H_ei_mean(1) = start.H_ei_mean;
    H_ie_mean(1) = start.H_ie_mean;
    H_ii_mean(1) = start.H_ii_mean;
    H_ee_var(1) = start.H_ee_var;
    H_ei_var(1) = start.H_ei_var;
    H_ie_var(1) = start.H_ie_var;
    H_ii_var(1) = start.H_ii_var;
    
    n_e(1,:) = start.n_e;
    n_i(1,:) = start.n_i;
    V_e_all = start.V_e_all;
    V_i_all = start.V_i_all;
    V_e_mean(1,:) = V_e_all./n_e(1,:);
    V_i_mean(1,:) = V_i_all./n_i(1,:);
    fr_e(1) = start.fr_e(1);
    fr_i(1) = start.fr_i(1);
    ref_e(1) = start.ref_e;
    ref_i(1) = start.ref_i;
end

for i = 2:params.t_end % 从“这一时刻”计算
    %% 将上一时刻的各种参数输入module中计算变化
    E_state.V_n_all = V_e_all(i-1,:);
    E_state.n_n = n_e(i-1,:);
    E_state.fr_n = fr_e(i-1);
    E_state.ref_n = ref_e(i-1);    
    I_ee_mean = params.s_ee*H_ee_mean(i-1)/params.tau_ee;
    I_ei_mean = params.s_ei*H_ei_mean(i-1)/params.tau_i*(V_e_mean(i-1,:)+params.Mr)/(params.M+params.Mr);
    I_ee_var = params.s_ee*H_ee_var(i-1)/params.tau_ee;
    I_ei_var = params.s_ei*H_ei_var(i-1)/params.tau_i*(V_e_mean(i-1,:)+params.Mr)/(params.M+params.Mr);
    I_e_mean(i,:) = I_eex(i-1,:) + I_ee_mean - I_ei_mean;
    I_e_var(i,:) = I_eex(i-1,:) + I_ee_var + I_ei_var;
    E_state.I_n_mean = I_e_mean(i,:);
    E_state.I_n_var = I_e_var(i,:);
    E_output = DIFODE_module(E_state,params);
    
    I_state.V_n_all = V_i_all(i-1,:);
    I_state.n_n = n_i(i-1,:);
    I_state.fr_n = fr_i(i-1);
    I_state.ref_n = ref_i(i-1);    
    I_ie_mean = params.s_ie*H_ie_mean(i-1)/params.tau_ie;
    I_ii_mean = params.s_ii*H_ii_mean(i-1)/params.tau_i*(V_i_mean(i-1,:)+params.Mr)/(params.M+params.Mr);
    I_ie_var = params.s_ie*H_ie_var(i-1)/params.tau_ie;
    I_ii_var = params.s_ii*H_ei_var(i-1)/params.tau_i*(V_i_mean(i-1,:)+params.Mr)/(params.M+params.Mr);
    I_i_mean(i,:) = I_iex(i-1,:) + I_ie_mean - I_ii_mean;
    I_i_var(i,:) = I_iex(i-1,:) + I_ie_var + I_ii_var;
    I_state.I_n_mean = I_i_mean(i,:);
    I_state.I_n_var = I_i_var(i,:);
    I_output = DIFODE_module(I_state,params);
    
    %% 计算这一时刻电流的分布参数
    dH_ee_mean = -H_ee_mean(i-1)/params.tau_ee + fr_e(i-1)*params.p_ee; % mean
    dH_ee_var = -H_ee_var(i-1)*2/(params.tau_ee) + H_ee_mean(i-1)/params.tau_ee + fr_e(i-1)*params.p_ee*(1-params.p_ee); % var
    H_ee_mean(i) = H_ee_mean(i-1) + dH_ee_mean*params.dt;
    H_ee_var(i) = H_ee_var(i-1) + dH_ee_var*params.dt;
    dH_ei_mean = -H_ei_mean(i-1)/params.tau_i + fr_i(i-1)*params.p_ei; % mean
    dH_ei_var = -H_ei_var(i-1)*2/(params.tau_i) + H_ei_mean(i-1)/params.tau_i + fr_i(i-1)*params.p_ei*(1-params.p_ei); % var
    H_ei_mean(i) = H_ei_mean(i-1) + dH_ei_mean*params.dt;
    H_ei_var(i) = H_ei_var(i-1) + dH_ei_var*params.dt;
    
    dH_ie_mean = -H_ie_mean(i-1)/params.tau_ie + fr_e(i-1)*params.p_ie; % mean
    dH_ie_var = -H_ie_var(i-1)*2/(params.tau_ie) + H_ie_mean(i-1)/params.tau_ie + fr_e(i-1)*params.p_ie*(1-params.p_ie); % var
    H_ie_mean(i) = H_ie_mean(i-1) + dH_ie_mean*params.dt;
    H_ie_var(i) = H_ie_var(i-1) + dH_ie_var*params.dt;
    dH_ii_mean = -H_ii_mean(i-1)/params.tau_i + fr_i(i-1)*params.p_ii; % mean
    dH_ii_var = -H_ii_var(i-1)*2/(params.tau_i) + H_ii_mean(i-1)/params.tau_i + fr_i(i-1)*params.p_ii*(1-params.p_ii); % var
    H_ii_mean(i) = H_ii_mean(i-1) + dH_ii_mean*params.dt;
    H_ii_var(i) = H_ii_var(i-1) + dH_ii_var*params.dt;
    
    %% 将output变为input
    V_e_all(i,:) = E_output.V_n_all;
    V_e_mean(i,:) = E_output.V_n_mean;
    n_e(i,:) = E_output.n_n;
    fr_e(i) = E_output.fr_n;
    ref_e(i) = E_output.ref_n;
    
    V_i_all(i,:) = I_output.V_n_all;
    V_i_mean(i,:) = I_output.V_n_mean;
    n_i(i,:) = I_output.n_n;
    fr_i(i) = I_output.fr_n;
    ref_i(i) = I_output.ref_n;
end

res_DIFODE.n_e = n_e;
res_DIFODE.V_e_mean = V_e_mean;
res_DIFODE.V_e_all = V_e_all;
res_DIFODE.ref_e = ref_e;
res_DIFODE.H_ee_mean = H_ee_mean;
res_DIFODE.H_ei_mean = H_ei_mean;
res_DIFODE.H_ee_var = H_ee_var;
res_DIFODE.H_ei_var = H_ei_var;
res_DIFODE.fr_e = fr_e;
res_DIFODE.I_eex = I_eex;
res_DIFODE.I_e_mean = I_e_mean;
res_DIFODE.I_e_var = I_e_var;

res_DIFODE.fr_i = fr_i;
res_DIFODE.n_i = n_i;
res_DIFODE.V_i_mean = V_i_mean;
res_DIFODE.V_i_all = V_i_all;
res_DIFODE.ref_i = ref_i;
res_DIFODE.H_ie_mean = H_ie_mean;
res_DIFODE.H_ii_mean = H_ii_mean;
res_DIFODE.H_ie_var = H_ie_var;
res_DIFODE.H_ii_var = H_ii_var;
res_DIFODE.I_iex = I_iex;
res_DIFODE.I_i_mean = I_i_mean;
res_DIFODE.I_i_var = I_i_var;

clear n_e
clear V_e_all
clear V_e_mean
clear ref_e
clear H_ee_mean
clear H_ei_mean
clear H_ee_var
clear H_ei_var
clear dH_ee_mean
clear dH_ei_mean
clear dH_ee_var
clear dH_ei_var
clear I_ee_mean
clear I_ei_mean
clear I_ee_var
clear I_ei_var
clear I_e_mean
clear fr_e
clear I_e_var
clear E_output
clear I_eex

clear n_i
clear V_i_all
clear V_i_mean
clear ref_i
clear H_ie_mean
clear H_ii_mean
clear H_ie_var
clear H_ii_var
clear dH_ie_mean
clear dH_ii_mean
clear dH_ie_var
clear dH_ii_var
clear I_ie_mean
clear I_ii_mean
clear I_ie_var
clear I_ii_var
clear I_i_mean
clear fr_i
clear I_i_var
clear I_output
clear I_iex

clear i
toc





