function [J_mat,D_vac] = get_jacobian_matrix(res_ODE,params,time_index)
%% params
dt = 0.05;
dV = 0.01;
dx = 0.01;
dxn = 0.1;
ne = params.ne;
ni = params.ni;
pee = params.p_ee; 
pie = params.p_ie; 
pei = params.p_ei; 
pii = params.p_ii; 
tau_ee = params.tau_ee;
tau_ie = params.tau_ie;
tau_i = params.tau_i;
tau_r = params.tau_r;
V_bin = 10;
V_bin_min = -3;
V_bin_num = params.M/V_bin - V_bin_min;


%% 结果
start.H_ee_mean = res_ODE.H_ee_mean(time_index);
start.H_ei_mean = res_ODE.H_ei_mean(time_index);
start.H_ie_mean = res_ODE.H_ie_mean(time_index);
start.H_ii_mean = res_ODE.H_ii_mean(time_index);
start.H_ee_var = res_ODE.H_ee_var(time_index);
start.H_ei_var = res_ODE.H_ei_var(time_index);
start.H_ie_var = res_ODE.H_ie_var(time_index);
start.H_ii_var = res_ODE.H_ii_var(time_index);
start.nfe = res_ODE.nfe(time_index);
start.nfi = res_ODE.nfi(time_index);
start.fre = res_ODE.fre(time_index);
start.fri = res_ODE.fri(time_index);
start.ref_e = res_ODE.ref_e(time_index);
start.ref_i = res_ODE.ref_i(time_index);
start.n_e = res_ODE.n_e(time_index,:);
start.n_i = res_ODE.n_i(time_index,:);
start.V_e_all = res_ODE.V_e_all(time_index,:);
start.V_i_all = res_ODE.V_i_all(time_index,:);
params.duration_time = 2*dt;
res_ODEfull_modelx = run_MIFIODE2_model(params,start);

for i = 1
    J_mat = zeros(V_bin_num*2*2+10,V_bin_num*2*2+10);
    D_vac = zeros(V_bin_num*2*2+10,1);
    %% 前8个函数可以直接求出来
%     % 第1行
%     J_mat(1,1) = -1/tau_ee;
%     J_mat(1,9) = pee;
%     % 第2行
%     J_mat(2,2) = -1/tau_ei;
%     J_mat(2,10) = pei;
%     % 第3行
%     J_mat(3,3) = -1/tau_ie;
%     J_mat(3,9) = pie;
%     % 第4行
%     J_mat(4,4) = -1/tau_ii;
%     J_mat(4,10) = pii;
% 
%     % 第5行
%     J_mat(5,5) = -2/tau_ee;
%     J_mat(5,1) = -1/tau_ee;
%     J_mat(5,9) = pee*(1-pee);
%     % 第6行
%     J_mat(6,6) = -2/tau_ei;
%     J_mat(6,2) = -1/tau_ei;
%     J_mat(6,10) = pei*(1-pei);
%     % 第7行
%     J_mat(7,7) = -2/tau_ie;
%     J_mat(7,3) = -1/tau_ie;
%     J_mat(7,9) = pie*(1-pie);
%     % 第8行
%     J_mat(8,8) = -2/tau_ii;
%     J_mat(8,4) = -1/tau_ii;
%     J_mat(8,10) = pii*(1-pii);

    % 第1行
    startdx = start;
    startdx.H_ee_mean = start.H_ee_mean + dx;
    res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
    D_vac(1) = (res_ODEfull_modelx.H_ee_mean(2)-res_ODEfull_modelx.H_ee_mean(1))/dt;
    J_mat(1,1) = ((res_ODEfull_modelxdx.H_ee_mean(2)-res_ODEfull_modelxdx.H_ee_mean(1)) - ...
        (res_ODEfull_modelx.H_ee_mean(2)-res_ODEfull_modelx.H_ee_mean(1)))/dt/dx;
    
    startdx = start;
    startdx.fre = start.fre + dxn;
    res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
    J_mat(1,9) = ((res_ODEfull_modelxdx.H_ee_mean(2)-res_ODEfull_modelxdx.H_ee_mean(1)) - ...
        (res_ODEfull_modelx.H_ee_mean(2)-res_ODEfull_modelx.H_ee_mean(1)))/dxn/dt;
    
    % 第2行
    startdx = start;
    startdx.H_ei_mean = start.H_ei_mean + dx;
    res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
    D_vac(2) = (res_ODEfull_modelx.H_ei_mean(2)-res_ODEfull_modelx.H_ei_mean(1))/dt;
    J_mat(2,2) = ((res_ODEfull_modelxdx.H_ei_mean(2)-res_ODEfull_modelxdx.H_ei_mean(1)) - ...
        (res_ODEfull_modelx.H_ei_mean(2)-res_ODEfull_modelx.H_ei_mean(1)))/dt/dx;
    
    startdx = start;
    startdx.fri = start.fri + dxn;
    res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
    J_mat(2,10) = ((res_ODEfull_modelxdx.H_ei_mean(2)-res_ODEfull_modelxdx.H_ei_mean(1)) - ...
        (res_ODEfull_modelx.H_ei_mean(2)-res_ODEfull_modelx.H_ei_mean(1)))/dt/dxn;
    
    % 第3行
    J_mat(3,3) = -1/tau_ie;
    J_mat(3,9) = pie;
    D_vac(3) = -1/tau_ie*res_ODEfull_modelx.H_ie_mean(1) + res_ODEfull_modelx.fre(1)*pie;
    
    % 第4行
    J_mat(4,4) = -1/tau_i;
    J_mat(4,10) = pii;
    D_vac(4) = -1/tau_i*res_ODEfull_modelx.H_ii_mean(1) + res_ODEfull_modelx.fri(1)*pii;

    % 第5行
    startdx = start;
    startdx.H_ee_mean = start.H_ee_mean + dx;
    res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
    D_vac(5) = res_ODEfull_modelx.H_ee_var(2)-res_ODEfull_modelx.H_ee_var(1);
    J_mat(5,1) = ((res_ODEfull_modelxdx.H_ee_var(2)-res_ODEfull_modelxdx.H_ee_var(1)) - ...
        (res_ODEfull_modelx.H_ee_var(2)-res_ODEfull_modelx.H_ee_var(1)))/dt/dx;
    
    startdx = start;
    startdx.H_ee_var = start.H_ee_var + dx;
    res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
    J_mat(5,5) = ((res_ODEfull_modelxdx.H_ee_var(2)-res_ODEfull_modelxdx.H_ee_var(1)) - ...
        (res_ODEfull_modelx.H_ee_var(2)-res_ODEfull_modelx.H_ee_var(1)))/dt/dx;
    
    startdx = start;
    startdx.fre = start.fre + dxn;
    res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
    J_mat(5,9) = ((res_ODEfull_modelxdx.H_ee_var(2)-res_ODEfull_modelxdx.H_ee_var(1)) - ...
        (res_ODEfull_modelx.H_ee_var(2)-res_ODEfull_modelx.H_ee_var(1)))/dt/dxn;
    
    % 第6行
    J_mat(6,6) = -2/tau_i;
    J_mat(6,2) = 1/tau_i;
    J_mat(6,10) = pei*(1-pei);
    D_vac(6) = -2/tau_i*res_ODEfull_modelx.H_ei_var(1) + 1/tau_i*res_ODEfull_modelx.H_ei_mean(1) + res_ODEfull_modelx.fri(1)*pei*(1-pei);
    
    % 第7行
    J_mat(7,7) = -2/tau_ie;
    J_mat(7,3) = 1/tau_ie;
    J_mat(7,9) = pie*(1-pie);
    D_vac(7) = -2/tau_ie*res_ODEfull_modelx.H_ie_var(1) + 1/tau_ie*res_ODEfull_modelx.H_ie_mean(1) + res_ODEfull_modelx.fre(1)*pie*(1-pei);
    
    % 第8行
    J_mat(8,8) = -2/tau_i;
    J_mat(8,4) = 1/tau_i;
    J_mat(8,10) = pii*(1-pii);
    D_vac(8) = -2/tau_i*res_ODEfull_modelx.H_ii_var(1) + 1/tau_i*res_ODEfull_modelx.H_ii_mean(1) + res_ODEfull_modelx.fri(1)*pii*(1-pii);
    
    for j = 1:V_bin_num
       %% E神经元H变量
        startdx = start;
        startdx.H_ee_mean = start.H_ee_mean + dx;
        res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
        D_vac(10+2*j-1) = (res_ODEfull_modelx.n_e(2,j)-res_ODEfull_modelx.n_e(1,j))/dt;
        J_mat(10+2*j-1,1) = ((res_ODEfull_modelxdx.n_e(2,j)-res_ODEfull_modelxdx.n_e(1,j)) - ...
            (res_ODEfull_modelx.n_e(2,j)-res_ODEfull_modelx.n_e(1,j)))/dt/dx;
        D_vac(10+2*j) = (res_ODEfull_modelx.V_e_all(2,j)-res_ODEfull_modelx.V_e_all(1,j))/dt;
        J_mat(10+2*j,1) = ((res_ODEfull_modelxdx.V_e_all(2,j)-res_ODEfull_modelxdx.V_e_all(1,j)) - ...
            (res_ODEfull_modelx.V_e_all(2,j)-res_ODEfull_modelx.V_e_all(1,j)))/dt/dx;

        startdx = start;
        startdx.H_ei_mean = start.H_ei_mean + dx;
        res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
        J_mat(10+2*j-1,2) = ((res_ODEfull_modelxdx.n_e(2,j)-res_ODEfull_modelxdx.n_e(1,j)) - ...
            (res_ODEfull_modelx.n_e(2,j)-res_ODEfull_modelx.n_e(1,j)))/dt/dx;
        J_mat(10+2*j,2) = ((res_ODEfull_modelxdx.V_e_all(2,j)-res_ODEfull_modelxdx.V_e_all(1,j)) - ...
            (res_ODEfull_modelx.V_e_all(2,j)-res_ODEfull_modelx.V_e_all(1,j)))/dt/dx;

        startdx = start;
        startdx.H_ee_var = start.H_ee_var + dx;
        res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
        J_mat(10+2*j-1,5) = ((res_ODEfull_modelxdx.n_e(2,j)-res_ODEfull_modelxdx.n_e(1,j)) - ...
            (res_ODEfull_modelx.n_e(2,j)-res_ODEfull_modelx.n_e(1,j)))/dt/dx;
        J_mat(10+2*j,5) = ((res_ODEfull_modelxdx.V_e_all(2,j)-res_ODEfull_modelxdx.V_e_all(1,j)) - ...
            (res_ODEfull_modelx.V_e_all(2,j)-res_ODEfull_modelx.V_e_all(1,j)))/dt/dx;

        startdx = start;
        startdx.H_ei_var = start.H_ei_var + dx;
        res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
        J_mat(10+2*j-1,6) = ((res_ODEfull_modelxdx.n_e(2,j)-res_ODEfull_modelxdx.n_e(1,j)) - ...
            (res_ODEfull_modelx.n_e(2,j)-res_ODEfull_modelx.n_e(1,j)))/dt/dx;
        J_mat(10+2*j,6) = ((res_ODEfull_modelxdx.V_e_all(2,j)-res_ODEfull_modelxdx.V_e_all(1,j)) - ...
            (res_ODEfull_modelx.V_e_all(2,j)-res_ODEfull_modelx.V_e_all(1,j)))/dt/dx;

        % 自己的区间
        if start.n_e(j) ~= 0
            startdx = start;
            startdx.n_e(j) = start.n_e(j) + dxn;% 这里的n不是rate，所以不用/dt
            res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
            J_mat(10+2*j-1,10+2*j-1) = ((res_ODEfull_modelxdx.n_e(2,j)-res_ODEfull_modelxdx.n_e(1,j)) - ...
                (res_ODEfull_modelx.n_e(2,j)-res_ODEfull_modelx.n_e(1,j)))/dt/dxn;
            J_mat(10+2*j,10+2*j-1) = ((res_ODEfull_modelxdx.V_e_all(2,j)-res_ODEfull_modelxdx.V_e_all(1,j)) - ...
                (res_ODEfull_modelx.V_e_all(2,j)-res_ODEfull_modelx.V_e_all(1,j)))/dt/dxn;

            startdx = start;
            startdx.V_e_all(j) = start.V_e_all(j) + dx;
            res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
            J_mat(10+2*j-1,10+2*j) = ((res_ODEfull_modelxdx.n_e(2,j)-res_ODEfull_modelxdx.n_e(1,j)) - ...
                (res_ODEfull_modelx.n_e(2,j)-res_ODEfull_modelx.n_e(1,j)))/dx/dt;
            J_mat(10+2*j,10+2*j) = ((res_ODEfull_modelxdx.V_e_all(2,j)-res_ODEfull_modelxdx.V_e_all(1,j)) - ...
                (res_ODEfull_modelx.V_e_all(2,j)-res_ODEfull_modelx.V_e_all(1,j)))/dx/dt;
        end
        
        %% I神经元H变量
        startdx = start;
        startdx.H_ie_mean = start.H_ie_mean + dx;
        res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
        D_vac(10+2*V_bin_num+2*j-1) = (res_ODEfull_modelx.n_i(2,j)-res_ODEfull_modelx.n_i(1,j))/dt;
        J_mat(10+2*V_bin_num+2*j-1,3) = ((res_ODEfull_modelxdx.n_i(2,j)-res_ODEfull_modelxdx.n_i(1,j)) - ...
            (res_ODEfull_modelx.n_i(2,j)-res_ODEfull_modelx.n_i(1,j)))/dt/dx;
        D_vac(10+2*V_bin_num+2*j) = (res_ODEfull_modelx.V_i_all(2,j)-res_ODEfull_modelx.V_i_all(1,j))/dt;
        J_mat(10+2*V_bin_num+2*j,3) = ((res_ODEfull_modelxdx.V_i_all(2,j)-res_ODEfull_modelxdx.V_i_all(1,j)) - ...
            (res_ODEfull_modelx.V_i_all(2,j)-res_ODEfull_modelx.V_i_all(1,j)))/dt/dx;
        
        startdx = start;
        startdx.H_ii_mean = start.H_ii_mean+dx;
        res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
        J_mat(10+2*V_bin_num+2*j-1,4) = ((res_ODEfull_modelxdx.n_i(2,j)-res_ODEfull_modelxdx.n_i(1,j)) - ...
            (res_ODEfull_modelx.n_i(2,j)-res_ODEfull_modelx.n_i(1,j)))/dt/dx;
        J_mat(10+2*V_bin_num+2*j,4) = ((res_ODEfull_modelxdx.V_i_all(2,j)-res_ODEfull_modelxdx.V_i_all(1,j)) - ...
            (res_ODEfull_modelx.V_i_all(2,j)-res_ODEfull_modelx.V_i_all(1,j)))/dt/dx;

        startdx = start;
        startdx.H_ie_var = start.H_ie_var+dx;
        res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
        J_mat(10+2*V_bin_num+2*j-1,7) = ((res_ODEfull_modelxdx.n_i(2,j)-res_ODEfull_modelxdx.n_i(1,j)) - ...
            (res_ODEfull_modelx.n_i(2,j)-res_ODEfull_modelx.n_i(1,j)))/dt/dx;
        J_mat(10+2*V_bin_num+2*j,7) = ((res_ODEfull_modelxdx.V_i_all(2,j)-res_ODEfull_modelxdx.V_i_all(1,j)) - ...
            (res_ODEfull_modelx.V_i_all(2,j)-res_ODEfull_modelx.V_i_all(1,j)))/dt/dx;

        startdx = start;
        startdx.H_ii_var = start.H_ii_var+dx;
        res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
        J_mat(10+2*V_bin_num+2*j-1,8) = ((res_ODEfull_modelxdx.n_i(2,j)-res_ODEfull_modelxdx.n_i(1,j)) - ...
            (res_ODEfull_modelx.n_i(2,j)-res_ODEfull_modelx.n_i(1,j)))/dt/dx;
        J_mat(10+2*V_bin_num+2*j,8) = ((res_ODEfull_modelxdx.V_i_all(2,j)-res_ODEfull_modelxdx.V_i_all(1,j)) - ...
            (res_ODEfull_modelx.V_i_all(2,j)-res_ODEfull_modelx.V_i_all(1,j)))/dt/dx;

        % 自己的区间
        if start.n_i(j) ~= 0
            startdx = start;
            startdx.n_i(j) = start.n_i(j)+dxn;
            res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
            J_mat(10+2*V_bin_num+2*j-1,10+2*V_bin_num+2*j-1) = ((res_ODEfull_modelxdx.n_i(2,j)-res_ODEfull_modelxdx.n_i(1,j)) - ...
                (res_ODEfull_modelx.n_i(2,j)-res_ODEfull_modelx.n_i(1,j)))/dt/dxn;
            J_mat(10+2*V_bin_num+2*j,10+2*V_bin_num+2*j-1) = ((res_ODEfull_modelxdx.V_i_all(2,j)-res_ODEfull_modelxdx.V_i_all(1,j)) - ...
                (res_ODEfull_modelx.V_i_all(2,j)-res_ODEfull_modelx.V_i_all(1,j)))/dt/dxn;

            startdx = start;
            startdx.V_i_all(j) = start.V_i_all(j)+dx;
            res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
            J_mat(10+2*V_bin_num+2*j-1,10+2*V_bin_num+2*j) = ((res_ODEfull_modelxdx.n_i(2,j)-res_ODEfull_modelxdx.n_i(1,j)) - ...
                (res_ODEfull_modelx.n_i(2,j)-res_ODEfull_modelx.n_i(1,j)))/dt/dx;
            J_mat(10+2*V_bin_num+2*j,10+2*V_bin_num+2*j) = ((res_ODEfull_modelxdx.V_i_all(2,j)-res_ODEfull_modelxdx.V_i_all(1,j)) - ...
                (res_ODEfull_modelx.V_i_all(2,j)-res_ODEfull_modelx.V_i_all(1,j)))/dt/dx;
        end
        
        if j==1
            %% E神经元右侧落下的区间
            if start.n_e(j+1) ~= 0
                startdx = start;
                startdx.n_e(j+1) = start.n_e(j+1)+dxn;
                res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
                J_mat(10+2*j-1,10+2*j+1) = ((res_ODEfull_modelxdx.n_e(2,j)-res_ODEfull_modelxdx.n_e(1,j)) - ...
                    (res_ODEfull_modelx.n_e(2,j)-res_ODEfull_modelx.n_e(1,j)))/dt/dxn;
                J_mat(10+2*j,10+2*j+1) = ((res_ODEfull_modelxdx.V_e_all(2,j)-res_ODEfull_modelxdx.V_e_all(1,j)) - ...
                    (res_ODEfull_modelx.V_e_all(2,j)-res_ODEfull_modelx.V_e_all(1,j)))/dt/dxn;

                startdx = start;
                startdx.V_e_all(j+1) = start.V_e_all(j+1)+dx;
                res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
                J_mat(10+2*j-1,10+2*j+2) = ((res_ODEfull_modelxdx.n_e(2,j)-res_ODEfull_modelxdx.n_e(1,j)) - ...
                    (res_ODEfull_modelx.n_e(2,j)-res_ODEfull_modelx.n_e(1,j)))/dt/dx;
                J_mat(10+2*j,10+2*j+2) = ((res_ODEfull_modelxdx.V_e_all(2,j)-res_ODEfull_modelxdx.V_e_all(1,j)) - ...
                    (res_ODEfull_modelx.V_e_all(2,j)-res_ODEfull_modelx.V_e_all(1,j)))/dt/dx;
            end
            
            %% I神经元右侧落下的区间
            if start.n_i(j+1) ~= 0
                startdx = start;
                startdx.n_i(j+1) = start.n_i(j+1)+dxn;
                res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
                J_mat(10+2*V_bin_num+2*j-1,10+2*V_bin_num+2*j+1) = ((res_ODEfull_modelxdx.n_i(2,j)-res_ODEfull_modelxdx.n_i(2,j)) - ...
                    (res_ODEfull_modelx.n_i(2,j)-res_ODEfull_modelx.n_i(2,j)))/dt/dxn;
                J_mat(10+2*V_bin_num+2*j,10+2*V_bin_num+2*j+1) = ((res_ODEfull_modelxdx.V_i_all(2,j)-res_ODEfull_modelxdx.V_i_all(2,j)) - ...
                    (res_ODEfull_modelx.V_i_all(2,j)-res_ODEfull_modelx.V_i_all(2,j)))/dt/dxn;

                startdx = start;
                startdx.V_i_all(j+1) = start.V_i_all(j+1)+dx;
                res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
                J_mat(10+2*V_bin_num+2*j-1,10+2*V_bin_num+2*j+2) = ((res_ODEfull_modelxdx.n_i(2,j)-res_ODEfull_modelxdx.n_i(2,j)) - ...
                    (res_ODEfull_modelx.n_i(2,j)-res_ODEfull_modelx.n_i(2,j)))/dt/dx;
                J_mat(10+2*V_bin_num+2*j,10+2*V_bin_num+2*j+2) = ((res_ODEfull_modelxdx.V_i_all(2,j)-res_ODEfull_modelxdx.V_i_all(2,j)) - ...
                    (res_ODEfull_modelx.V_i_all(2,j)-res_ODEfull_modelx.V_i_all(2,j)))/dt/dx;
            end
            
        elseif j==V_bin_num
            %% E神经元fre
            startdx = start;
            startdx.H_ee_mean = start.H_ee_mean + dx;
            res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
            D_vac(9) = (res_ODEfull_modelx.fre(2) - res_ODEfull_modelx.fre(1))/dt;
            J_mat(9,1) = ((res_ODEfull_modelxdx.fre(2) - res_ODEfull_modelxdx.fre(1)) - ...
                (res_ODEfull_modelx.fre(2) - res_ODEfull_modelx.fre(1)))/dt/dx;

            startdx = start;
            startdx.H_ei_mean = start.H_ei_mean + dx;
            res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
            J_mat(9,2) = ((res_ODEfull_modelxdx.fre(2) - res_ODEfull_modelxdx.fre(1)) - ...
                (res_ODEfull_modelx.fre(2) - res_ODEfull_modelx.fre(1)))/dt/dx;

            startdx = start;
            startdx.H_ee_var = start.H_ee_var + dx;
            res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
            J_mat(9,5) = ((res_ODEfull_modelxdx.fre(2) - res_ODEfull_modelxdx.fre(1)) - ...
                (res_ODEfull_modelx.fre(2) - res_ODEfull_modelx.fre(1)))/dt/dx;

            startdx = start;
            startdx.H_ei_var = start.H_ei_var + dx;
            res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
            J_mat(9,6) = ((res_ODEfull_modelxdx.fre(2) - res_ODEfull_modelxdx.fre(1)) - ...
                (res_ODEfull_modelx.fre(2) - res_ODEfull_modelx.fre(1)))/dt/dx;
            
            if start.n_e(j) ~= 0
                startdx = start;
                startdx.n_e(j) = start.n_e(j)+dxn;
                res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
                J_mat(9,10+2*j-1) = ((res_ODEfull_modelxdx.fre(2) - res_ODEfull_modelxdx.fre(1)) - ...
                    (res_ODEfull_modelx.fre(2) - res_ODEfull_modelx.fre(1)))/dt/dxn;
                
                startdx = start;
                startdx.V_e_all(j) = start.V_e_all(j)+dx;
                res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
                J_mat(9,10+2*j) = ((res_ODEfull_modelxdx.fre(2) - res_ODEfull_modelxdx.fre(1)) - ...
                    (res_ODEfull_modelx.fre(2) - res_ODEfull_modelx.fre(1)))/dt/dx;
            end
            
            %% I神经元fri
            startdx = start;
            startdx.H_ie_mean = start.H_ie_mean + dx;
            res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
            D_vac(10) = (res_ODEfull_modelx.fri(2) - res_ODEfull_modelx.fri(1))/dt;
            J_mat(10,3) = ((res_ODEfull_modelxdx.fri(2) - res_ODEfull_modelxdx.fri(1)) - ...
                (res_ODEfull_modelx.fri(2) - res_ODEfull_modelx.fri(1)))/dt/dx;

            startdx = start;
            startdx.H_ii_mean = start.H_ii_mean + dx;
            res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
            J_mat(10,4) = ((res_ODEfull_modelxdx.fri(2) - res_ODEfull_modelxdx.fri(1)) - ...
                (res_ODEfull_modelx.fri(2) - res_ODEfull_modelx.fri(1)))/dt/dx;

            startdx = start;
            startdx.H_ie_var = start.H_ie_var + dx;
            res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
            J_mat(10,7) = ((res_ODEfull_modelxdx.fri(2) - res_ODEfull_modelxdx.fri(1)) - ...
                (res_ODEfull_modelx.fri(2) - res_ODEfull_modelx.fri(1)))/dt/dx;

            startdx = start;
            startdx.H_ii_var = start.H_ii_var + dx;
            res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
            J_mat(10,8) = ((res_ODEfull_modelxdx.fri(2) - res_ODEfull_modelxdx.fri(1)) - ...
                (res_ODEfull_modelx.fri(2) - res_ODEfull_modelx.fri(1)))/dt/dx;
            
            if start.n_i(j) ~= 0
                startdx = start;
                startdx.n_i(j) = start.n_i(j)+dxn;
                res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
                J_mat(10,10+2*V_bin_num+2*j-1) = ((res_ODEfull_modelxdx.fri(2) - res_ODEfull_modelxdx.fri(1)) - ...
                    (res_ODEfull_modelx.fri(2) - res_ODEfull_modelx.fri(1)))/dt/dxn;
                
                startdx = start;
                startdx.V_i_all(j) = start.V_i_all(j)+dx;
                res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
                J_mat(10,10+2*V_bin_num+2*j) = ((res_ODEfull_modelxdx.fri(2) - res_ODEfull_modelxdx.fri(1)) - ...
                    (res_ODEfull_modelx.fri(2) - res_ODEfull_modelx.fri(1)))/dt/dx;
            end
            
            %% E神经元左侧升上来的区间
            if start.n_e(j-1) ~= 0
                startdx = start;
                startdx.n_e(j-1) = start.n_e(j-1)+dxn;
                res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
                J_mat(10+2*j-1,10+2*j-3) = ((res_ODEfull_modelxdx.n_e(2,j)-res_ODEfull_modelxdx.n_e(1,j)) - ...
                    (res_ODEfull_modelx.n_e(2,j)-res_ODEfull_modelx.n_e(1,j)))/dt/dxn;
                J_mat(10+2*j,10+2*j-3) = ((res_ODEfull_modelxdx.V_e_all(2,j)-res_ODEfull_modelxdx.V_e_all(1,j)) - ...
                    (res_ODEfull_modelx.V_e_all(2,j)-res_ODEfull_modelx.V_e_all(1,j)))/dt/dxn;
            
                startdx = start;
                startdx.V_e_all(j-1) = start.V_e_all(j-1)+dx;
                res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
                J_mat(10+2*j-1,10+2*j-2) = ((res_ODEfull_modelxdx.n_e(2,j)-res_ODEfull_modelxdx.n_e(1,j)) - ...
                    (res_ODEfull_modelx.n_e(2,j)-res_ODEfull_modelx.n_e(1,j)))/dt/dx;
                J_mat(10+2*j,10+2*j-2) = ((res_ODEfull_modelxdx.V_e_all(2,j)-res_ODEfull_modelxdx.V_e_all(1,j)) - ...
                    (res_ODEfull_modelx.V_e_all(2,j)-res_ODEfull_modelx.V_e_all(1,j)))/dt/dx;
            end
            
            %% I神经元左侧升上来的区间
            if start.n_i(j-1) ~= 0
                startdx = start;
                startdx.n_i(j-1) = start.n_i(j-1)+dxn;
                res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
                J_mat(10+2*V_bin_num+2*j-1,10+2*V_bin_num+2*j-3) = ((res_ODEfull_modelxdx.n_i(2,j)-res_ODEfull_modelxdx.n_i(1,j)) - ...
                    (res_ODEfull_modelx.n_i(2,j)-res_ODEfull_modelx.n_i(1,j)))/dt/dxn;
                J_mat(10+2*V_bin_num+2*j,10+2*V_bin_num+2*j-3) = ((res_ODEfull_modelxdx.V_i_all(2,j)-res_ODEfull_modelxdx.V_i_all(1,j)) - ...
                    (res_ODEfull_modelx.V_i_all(2,j)-res_ODEfull_modelx.V_i_all(1,j)))/dt/dxn;

                startdx = start;
                startdx.V_i_all(j-1) = start.V_i_all(j-1)+dx;
                res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
                J_mat(10+2*V_bin_num+2*j-1,10+2*V_bin_num+2*j-2) = ((res_ODEfull_modelxdx.n_i(2,j)-res_ODEfull_modelxdx.n_i(1,j)) - ...
                    (res_ODEfull_modelx.n_i(2,j)-res_ODEfull_modelx.n_i(1,j)))/dt/dx;
                J_mat(10+2*V_bin_num+2*j,10+2*V_bin_num+2*j-2) = ((res_ODEfull_modelxdx.V_i_all(2,j)-res_ODEfull_modelxdx.V_i_all(1,j)) - ...
                    (res_ODEfull_modelx.V_i_all(2,j)-res_ODEfull_modelx.V_i_all(1,j)))/dt/dx;
            end
            
        else % 其他所有区间  
            %% E神经元左侧升上来的区间
            if start.n_e(j-1) ~= 0
                startdx = start;
                startdx.n_e(j-1) = start.n_e(j-1)+dxn;
                res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
                J_mat(10+2*j-1,10+2*j-3) = ((res_ODEfull_modelxdx.n_e(2,j)-res_ODEfull_modelxdx.n_e(1,j)) - ...
                    (res_ODEfull_modelx.n_e(2,j)-res_ODEfull_modelx.n_e(1,j)))/dt/dxn;
                J_mat(10+2*j,10+2*j-3) = ((res_ODEfull_modelxdx.V_e_all(2,j)-res_ODEfull_modelxdx.V_e_all(1,j)) - ...
                    (res_ODEfull_modelx.V_e_all(2,j)-res_ODEfull_modelx.V_e_all(1,j)))/dt/dxn;

                startdx = start;
                startdx.V_e_all(j-1) = start.V_e_all(j-1)+dx;
                res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
                J_mat(10+2*j-1,10+2*j-2) = ((res_ODEfull_modelxdx.n_e(2,j)-res_ODEfull_modelxdx.n_e(1,j)) - ...
                    (res_ODEfull_modelx.n_e(2,j)-res_ODEfull_modelx.n_e(1,j)))/dt/dx;
                J_mat(10+2*j,10+2*j-2) = ((res_ODEfull_modelxdx.V_e_all(2,j)-res_ODEfull_modelxdx.V_e_all(1,j)) - ...
                    (res_ODEfull_modelx.V_e_all(2,j)-res_ODEfull_modelx.V_e_all(1,j)))/dt/dx;
            end

            
            %% E神经元右侧落上来的区间
            if start.n_e(j+1) ~= 0
                startdx = start;
                startdx.n_e(j+1) = start.n_e(j+1)+dxn;
                res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
                J_mat(10+2*j-1,10+2*j+1) = ((res_ODEfull_modelxdx.n_e(2,j)-res_ODEfull_modelxdx.n_e(1,j)) - ...
                    (res_ODEfull_modelx.n_e(2,j)-res_ODEfull_modelx.n_e(1,j)))/dt/dxn;
                J_mat(10+2*j,10+2*j+1) = ((res_ODEfull_modelxdx.V_e_all(2,j)-res_ODEfull_modelxdx.V_e_all(1,j)) - ...
                    (res_ODEfull_modelx.V_e_all(2,j)-res_ODEfull_modelx.V_e_all(1,j)))/dt/dxn;

                startdx = start;
                startdx.V_e_all(j+1) = start.V_e_all(j+1)+dx;
                res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
                J_mat(10+2*j-1,10+2*j+2) = ((res_ODEfull_modelxdx.n_e(2,j)-res_ODEfull_modelxdx.n_e(1,j)) - ...
                    (res_ODEfull_modelx.n_e(2,j)-res_ODEfull_modelx.n_e(1,j)))/dt/dx;
                J_mat(10+2*j,10+2*j+2) = ((res_ODEfull_modelxdx.V_e_all(2,j)-res_ODEfull_modelxdx.V_e_all(1,j)) - ...
                    (res_ODEfull_modelx.V_e_all(2,j)-res_ODEfull_modelx.V_e_all(1,j)))/dt/dx;
            end
            
            %% I神经元左侧升上来的区间
            if start.n_i(j-1) ~= 0
                startdx = start;
                startdx.n_i(j-1) = start.n_i(j-1)+dxn;
                res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
                J_mat(10+2*V_bin_num+2*j-1,10+2*V_bin_num+2*j-3) = ((res_ODEfull_modelxdx.n_i(2,j)-res_ODEfull_modelxdx.n_i(1,j)) - ...
                    (res_ODEfull_modelx.n_i(2,j)-res_ODEfull_modelx.n_i(1,j)))/dt/dxn;
                J_mat(10+2*V_bin_num+2*j,10+2*V_bin_num+2*j-3) = ((res_ODEfull_modelxdx.V_i_all(2,j)-res_ODEfull_modelxdx.V_i_all(1,j)) - ...
                    (res_ODEfull_modelx.V_i_all(2,j)-res_ODEfull_modelx.V_i_all(1,j)))/dt/dxn;

                startdx = start;
                startdx.V_i_all(j-1) = start.V_i_all(j-1)+dx;
                res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
                J_mat(10+2*V_bin_num+2*j-1,10+2*V_bin_num+2*j-2) = ((res_ODEfull_modelxdx.n_i(2,j)-res_ODEfull_modelxdx.n_i(1,j)) - ...
                    (res_ODEfull_modelx.n_i(2,j)-res_ODEfull_modelx.n_i(1,j)))/dt/dx;
                J_mat(10+2*V_bin_num+2*j,10+2*V_bin_num+2*j-2) = ((res_ODEfull_modelxdx.V_i_all(2,j)-res_ODEfull_modelxdx.V_i_all(1,j)) - ...
                    (res_ODEfull_modelx.V_i_all(2,j)-res_ODEfull_modelx.V_i_all(1,j)))/dt/dx;
            end
            
            %% I神经元右侧落下来的区间
            if start.n_i(j+1) ~= 0
                startdx = start;
                startdx.n_i(j+1) = start.n_i(j+1)+dxn;
                res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
                J_mat(10+2*V_bin_num+2*j-1,10+2*V_bin_num+2*j+1) = ((res_ODEfull_modelxdx.n_i(2,j)-res_ODEfull_modelxdx.n_i(1,j)) - ...
                    (res_ODEfull_modelx.n_i(2,j)-res_ODEfull_modelx.n_i(1,j)))/dt/dxn;
                J_mat(10+2*V_bin_num+2*j,10+2*V_bin_num+2*j+1) = ((res_ODEfull_modelxdx.V_i_all(2,j)-res_ODEfull_modelxdx.V_i_all(1,j)) - ...
                    (res_ODEfull_modelx.V_i_all(2,j)-res_ODEfull_modelx.V_i_all(1,j)))/dt/dxn;

                startdx = start;
                startdx.V_i_all(j+1) = start.V_i_all(j+1)+dx;
                res_ODEfull_modelxdx = run_MIFIODE2_model(params,startdx);
                J_mat(10+2*V_bin_num+2*j-1,10+2*V_bin_num+2*j+2) = ((res_ODEfull_modelxdx.n_i(2,j)-res_ODEfull_modelxdx.n_i(1,j)) - ...
                    (res_ODEfull_modelx.n_i(2,j)-res_ODEfull_modelx.n_i(1,j)))/dt/dx;
                J_mat(10+2*V_bin_num+2*j,10+2*V_bin_num+2*j+2) = ((res_ODEfull_modelxdx.V_i_all(2,j)-res_ODEfull_modelxdx.V_i_all(1,j)) - ...
                    (res_ODEfull_modelx.V_i_all(2,j)-res_ODEfull_modelx.V_i_all(1,j)))/dt/dx;
            end
        end
    end
    
    %% 标准化
    D_vac(1) = D_vac(1)/200;
    D_vac(2) = D_vac(2)/150;
    D_vac(3) = D_vac(3)/200;
    D_vac(4) = D_vac(4)/150;
    D_vac(5) = D_vac(1)/100;
    D_vac(6) = D_vac(2)/75;
    D_vac(7) = D_vac(3)/100;
    D_vac(8) = D_vac(4)/75;
    D_vac(9) = D_vac(9)/30;
    D_vac(10) = D_vac(10)/40;
    D_vac(11:2:(10+2*V_bin_num-1)) = D_vac(11:2:(10+2*V_bin_num-1))/120; 
    D_vac((10+2*V_bin_num+1):2:(10+4*V_bin_num-1)) = D_vac((10+2*V_bin_num+1):2:(10+4*V_bin_num-1))/40; 
    D_vac(12:2:(10+2*V_bin_num)) = D_vac(12:2:(10+2*V_bin_num))/params.M/120; 
    D_vac((10+2*V_bin_num+2):2:(10+4*V_bin_num)) = D_vac((10+2*V_bin_num+2):2:(10+4*V_bin_num))/params.M/40; 
    
    J_mat(1,:) = J_mat(1,:)/200;
    J_mat(2,:) = J_mat(2,:)/150;
    J_mat(3,:) = J_mat(3,:)/200;
    J_mat(4,:) = J_mat(4,:)/150;
    J_mat(5,:) = J_mat(5,:)/100;
    J_mat(6,:) = J_mat(6,:)/75;
    J_mat(7,:) = J_mat(7,:)/100;
    J_mat(8,:) = J_mat(8,:)/75;
    J_mat(9,:) = J_mat(9,:)/30;
    J_mat(10,:) = J_mat(10,:)/40;
    J_mat(11:2:(10+2*V_bin_num-1),:) = J_mat(11:2:(10+2*V_bin_num-1),:)/120;
    J_mat((10+2*V_bin_num+1):2:(10+4*V_bin_num-1),:) = J_mat((10+2*V_bin_num+1):2:(10+4*V_bin_num-1),:)/40; 
    J_mat(12:2:(10+2*V_bin_num),:) = J_mat(12:2:(10+2*V_bin_num),:)/params.M/120; 
    J_mat((10+2*V_bin_num+2):2:(10+4*V_bin_num),:) = J_mat((10+2*V_bin_num+2):2:(10+4*V_bin_num),:)/params.M/40; 

end
end

