function [IVrefMM_list] = generate_IVrefMM(estate_matref,istate_matref)
%% 参数
tau_ex = 1/7;
tau_ee = 1.4;
tau_ie = 1.2;
tau_i = 4.5;
tau_ref = 2;
ne=300;
ni=100;
pee=0.15;
pei=0.5;
pie=0.5;
pii=0.4;
s_ee = 4;
s_ie = 2;
s_ei = 4.91;
s_ii = 4.91;
M=100;
Mr=66;
duration_time = 1000;
min_V_e = -20;
max_V_e = 100;
min_V_i = -20;
max_V_i = 100;
min_I_e = -30;
max_I_e = 120;
min_I_i = -30;
max_I_i = 120;
refbin = 0.25;
max_ref_e = ceil(1/refbin);
min_ref_e = ceil(0/refbin);
max_ref_i = ceil(1/refbin);
min_ref_i = ceil(0/refbin);
V_e_statenum = max_V_e-min_V_e+1;
V_i_statenum = max_V_i-min_V_i+1;
I_e_statenum = max_I_e-min_I_e+1;
I_i_statenum = max_I_i-min_I_i+1;
ht = 0:0.1:10;
t_step = 0.1;


%% 分配空间
nfe = zeros(1,duration_time/t_step);
nfi = zeros(1,duration_time/t_step);
ne_ref = zeros(1,duration_time/t_step);
ni_ref = zeros(1,duration_time/t_step);
V_e_noref_mean = zeros(1,duration_time/t_step);
V_i_noref_mean = zeros(1,duration_time/t_step);
I_e_all_mean = zeros(1,duration_time/t_step);
I_i_all_mean = zeros(1,duration_time/t_step);


%% 初值
nfe(1) = 0;
nfi(1) = 0;
ne_ref(1) = nfe(1);
ni_ref(1) = nfi(1);
I_e_all_mean(1) = 7;
I_i_all_mean(1) = 7;
V_e_noref_mean(1) = 50;
V_i_noref_mean(1) = 40;


for i = 1:duration_time/t_step-1 % i从1开始计算，也就是i代表这一时刻
    %% 跳转下一刻的V mean变量与nf变量
     % E跳转
     V_e_index = V_e_noref_mean(i)-min_V_e+1; % 现在是要预测下一时刻的电压
     I_e_index = I_e_all_mean(i)-min_I_e+1;
     ne_ref_index = ceil(ne_ref(i)/ne/refbin)-min_ref_e+1;
     next_state_list = estate_matref{V_e_index,I_e_index,ne_ref_index};
     if isempty(next_state_list)
        next_state_list = vertcat(estate_matref{V_e_index-1:V_e_index+1,I_e_index-1:I_e_index+1,ne_ref_index});
     end
     if isempty(next_state_list)
        next_state_list = vertcat(estate_matref{V_e_index-3:V_e_index+3,I_e_index-3:I_e_index+3,max(1,ne_ref_index-1):min(5,ne_ref_index+1)});
     end
     if isempty(next_state_list)
        next_state_list = vertcat(estate_matref{V_e_index-8:V_e_index+8,I_e_index-8:I_e_index+8,max(1,ne_ref_index-2):min(5,ne_ref_index+2)});
     end
     state_index = randsrc(1,1,[1:length(next_state_list(:,3));next_state_list(:,3)'/sum(next_state_list(:,3))]);
     V_e_noref_mean(i+1) = next_state_list(state_index,1);
     nfe(i+1) = next_state_list(state_index,2);
     
     % I跳转
     V_i_index = V_i_noref_mean(i)-min_V_i+1; % 现在是要预测下一时刻的电压
     I_i_index = I_i_all_mean(i)-min_I_i+1;
     ni_ref_index = ceil(ni_ref(i)/ni/refbin)-min_ref_i+1;
     next_state_list = istate_matref{V_i_index,I_i_index,ni_ref_index};
     if isempty(next_state_list)
        next_state_list = vertcat(istate_matref{V_i_index-1:V_i_index+1,I_i_index-1:I_i_index+1,ni_ref_index});
     end
     if isempty(next_state_list)
        next_state_list = vertcat(istate_matref{V_i_index-3:V_i_index+3,I_i_index-3:I_i_index+3,max(1,ni_ref_index-1):min(5,ni_ref_index+1)});
     end
     if isempty(next_state_list)
        next_state_list = vertcat(istate_matref{V_i_index-8:V_i_index+8,I_i_index-8:I_i_index+8,max(1,ni_ref_index-2):min(5,ni_ref_index+2)});
     end
     state_index = randsrc(1,1,[1:length(next_state_list(:,3));next_state_list(:,3)'/sum(next_state_list(:,3))]);
     V_i_noref_mean(i+1) = next_state_list(state_index,1);
     nfi(i+1) = next_state_list(state_index,2);

     
    %% 计算下一刻的I
    h_tau_ee = exp(-(1/tau_ee)*ht);
    H_ee = pee * ne * conv(nfe,h_tau_ee);
    H_ee = H_ee(i+1);
    h_tau_ie = exp(-(1/tau_ie)*ht);
    H_ie = pie * ni * conv(nfe,h_tau_ie);
    H_ie = H_ie(i+1);

    h_tau_i = exp(-(1/tau_i)*ht);
    H_ei = pei * ne * conv(nfi,h_tau_i);
    H_ei = H_ei(i+1);
    h_tau_i = exp(-(1/tau_i)*ht);
    H_ii = pii * ni * conv(nfi,h_tau_i);
    H_ii = H_ii(i+1);
    
    h_tau_ref = exp(-(1/tau_ref)*ht);
    ref_e = conv(nfe,h_tau_ref);
    ne_ref(i+1) = round(ref_e(i+1));
    ref_i = conv(nfi,h_tau_ref);
    ni_ref(i+1) = round(ref_i(i+1));
    
    I_eex = 1/tau_ex*ne;
    I_iex = 1/tau_ex*ni;
    I_ee = s_ee*H_ee/tau_ee;
    I_ie = s_ie*H_ie/tau_ie;
    I_ei = (V_e_noref_mean(i+1)+Mr)*s_ei/(M+Mr)*H_ei/tau_i;
    I_ii = (V_i_noref_mean(i+1)+Mr)*s_ii/(M+Mr)*H_ii/tau_i;
    
    I_e_all_mean(i+1) = round((I_ee+I_eex-I_ei)/ne);
    I_i_all_mean(i+1) = round((I_ie+I_iex-I_ii)/ni);
     
     
end
    
    
IVrefMM_list = [nfe; nfi; V_e_noref_mean; V_i_noref_mean; I_e_all_mean; I_i_all_mean];
figure
plot(0.1:0.1:duration_time, IVrefMM_list(3,:),'m');
hold on
plot(0.1:0.1:duration_time, IVrefMM_list(4,:),'g');
legend('E V','I V');
xlabel('time(ms)');
ylabel('V');
title('IVMM V');

figure
plot(0.1:0.1:duration_time, IVrefMM_list(5,:),'r');
hold on
plot(0.1:0.1:duration_time, IVrefMM_list(6,:),'b');
legend('E current','I current');
xlabel('time(ms)');
ylabel('current');
title('IVMM current');

end


