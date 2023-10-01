function [res] = run_MIFIODE2_model(params,start)
res = struct;

%% preparation
dt = params.dt;
dV = 0.01;
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
digit_num = 2;

t_end = round(params.duration_time/dt);
V_e_mean = zeros(t_end, V_bin_num);
% V_e_mean(1,1-V_bin_min) = 5;
% V_e_mean(1,2-V_bin_min) = 15;
% V_e_mean(1,3-V_bin_min) = 25;
% V_e_mean(1,4-V_bin_min) = 35;
% V_e_mean(1,5-V_bin_min) = 45;

V_e_remain = zeros(t_end, V_bin_num);
V_e_up = zeros(t_end, V_bin_num);
V_e_down = zeros(t_end, V_bin_num);
I_e_mean = zeros(t_end, V_bin_num);
I_e_var = zeros(t_end, V_bin_num);

% V_i = repmat(-V_bin:V_bin:(M-V_bin),t_end,1);
V_i_mean = zeros(t_end, V_bin_num);
V_i_mean(1,1-V_bin_min) = 0;
V_i_remain = zeros(t_end, V_bin_num);
V_i_up = zeros(t_end, V_bin_num);
V_i_down = zeros(t_end, V_bin_num);
I_i_mean = zeros(t_end, V_bin_num);
I_i_var= zeros(t_end, V_bin_num);

n_e = zeros(t_end, V_bin_num);
n_e(1,1-V_bin_min) = ne;
% n_e(1,2-V_bin_min) = ne/5;
% n_e(1,3-V_bin_min) = ne/5;
% n_e(1,4-V_bin_min) = ne/5;
% n_e(1,5-V_bin_min) = ne/5;

n_e_remain = zeros(t_end, V_bin_num);
n_e_up = zeros(t_end, V_bin_num);
n_e_down = zeros(t_end, V_bin_num);
n_i = zeros(t_end, V_bin_num);
n_i(1,1-V_bin_min) = ni;
n_i_remain = zeros(t_end, V_bin_num);
n_i_up = zeros(t_end, V_bin_num+1);
n_i_down = zeros(t_end, V_bin_num+1);

nfe = zeros(t_end, 1);
nfi = zeros(t_end, 1);
fre = zeros(t_end, 1);
fri = zeros(t_end, 1);
ref_e = zeros(t_end, 1);
leave_e = zeros(t_end, 1);
leave_i = zeros(t_end, 1);
ref_i = zeros(t_end, 1);

H_ee_mean = zeros(t_end, 1);
H_ei_mean = zeros(t_end, 1);
H_ie_mean = zeros(t_end, 1);
H_ii_mean = zeros(t_end, 1);
H_ee_var = zeros(t_end, 1);
H_ei_var = zeros(t_end, 1);
H_ie_var = zeros(t_end, 1);
H_ii_var = zeros(t_end, 1);
if nargin == 2
    H_ee_mean(1) = start.H_ee_mean;
    H_ei_mean(1) = start.H_ei_mean;
    H_ie_mean(1) = start.H_ie_mean;
    H_ii_mean(1) = start.H_ii_mean;
    H_ee_var(1) = start.H_ee_var;
    H_ei_var(1) = start.H_ei_var;
    H_ie_var(1) = start.H_ie_var;
    H_ii_var(1) = start.H_ii_var;
    
    n_e(1,:) = start.n_e;
    n_e(1,end) = n_e(1,end);
    n_i(1,:) = start.n_i;
    V_e_all = start.V_e_all;
    V_i_all = start.V_i_all;
    V_e_mean(1,:) = V_e_all./n_e(1,:);
    V_i_mean(1,:) = V_i_all./n_i(1,:);
%     V_e_mean(1,:) = start.V_e_mean;
%     V_i_mean(1,:) = start.V_i_mean;
    fre(1) = start.fre(1);
    fri(1) = start.fri(1);
    nfe(1) = fre(1)*params.dt;
    nfi(1) = nfi(1)*params.dt;
%     nfi(1) = nfi(1)+6;

    ref_e(1) = start.ref_e;
    ref_e(1) = start.ref_e;    
    ref_i(1) = start.ref_i;
    
    for j = 1:V_bin_num    
        J_ee = params.s_ee;
        J_ei = (V_e_mean(1,j)+params.Mr)*params.s_ei/(params.M+params.Mr);    
        I_e_mean(1,j) = J_ee*H_ee_mean(1)/params.tau_ee - J_ei*H_ei_mean(1)/params.tau_i; % 注意这里有/dt有*dt，代表的含义不一样，所以分开写
        I_e_var(1,j) = J_ee*H_ee_var(1)/params.tau_ee + J_ei*H_ei_var(1)/params.tau_i;
        
        J_ie = params.s_ie;
        J_ii = (V_i_mean(1,j)+params.Mr)*params.s_ii/(params.M+params.Mr);    
        I_i_mean(1,j) = J_ie*H_ie_mean(1)/params.tau_ie - J_ii*H_ii_mean(1)/params.tau_i; % 注意这里有/dt有*dt，代表的含义不一样，所以分开写
        I_i_var(1,j) = J_ie*H_ie_var(1)/params.tau_ie + J_ii*H_ii_var(1)/params.tau_i;
    end
end
I_lowlim = -10;
I_highlim = 40;
IdeltaV = I_lowlim:dV:I_highlim;


%% run model
for i = 2:t_end % 从“这一时刻”计算
    nenow = round(sum(n_e(i-1,:))+ref_e(i-1));
    ninow = round(sum(n_i(i-1,:))+ref_i(i-1));
    if nenow~=300 || ninow~=100
        fprintf('error\n');
    end

    %% 计算每个区间神经元电流
    for j = 1:V_bin_num
         %% E神经元 
        if n_e(i-1,j) ~= 0 % 只有不为0的区间才有计算价值
            %% 计算V的概率分布（均匀分布）
            V_highlim = V_bin*(j+V_bin_min);
            V_lowlim = V_bin*(j+V_bin_min-1);
            if V_e_mean(i-1,j) < V_lowlim+dV
                V_e_p = zeros(1,(V_highlim-V_lowlim)/dV);
                V_e_p(1) = 1/dV;
            elseif V_e_mean(i-1,j) <= (V_highlim+V_lowlim)/2
                V_e_p = pdf('Uniform',(V_lowlim+dV):dV:V_highlim,V_lowlim,2*V_e_mean(i-1,j)-V_lowlim);
            elseif V_e_mean(i-1,j) >= V_highlim
                V_e_p = zeros(1,(V_highlim-V_lowlim)/dV);
                V_e_p(end) = 1/dV;
            else
                V_e_p = pdf('Uniform',(V_lowlim+dV):dV:V_highlim,2*V_e_mean(i-1,j)-V_highlim,V_highlim);
            end

            %% 计算这一时刻电流的分布参数以及上一时刻的电流分布（正态分布与泊松分布）
            dH_ee_mean = -H_ee_mean(i-1)/tau_ee + fre(i-1)*pee; % mean
            dH_ee_var = -H_ee_var(i-1)*2/(tau_ee) + H_ee_mean(i-1)/tau_ee + fre(i-1)*pee*(1-pee); % var
            H_ee_mean(i) = H_ee_mean(i-1) + dH_ee_mean*dt;
            H_ee_var(i) = H_ee_var(i-1) + dH_ee_var*dt;
            
            dH_ei_mean = -H_ei_mean(i-1)/tau_i + fri(i-1)*pei; % mean
            dH_ei_var = -H_ei_var(i-1)*2/(tau_i) + H_ei_mean(i-1)/tau_i + fri(i-1)*pei*(1-pei); % var
            H_ei_mean(i) = H_ei_mean(i-1) + dH_ei_mean*dt;
            H_ei_var(i) = H_ei_var(i-1) + dH_ei_var*dt;
            
            J_eex = 1;
            J_ee = params.s_ee;
            J_ei = (V_e_mean(i-1,j)+params.Mr)*params.s_ei/(params.M+params.Mr);    
            I_e_mean(i,j) = J_ee*H_ee_mean(i)/params.tau_ee - J_ei*H_ei_mean(i)/params.tau_i;
            I_e_var(i,j) = J_ee*H_ee_var(i)/params.tau_ee + J_ei*H_ei_var(i)/params.tau_i;
            I_e_p = pdf('Normal', IdeltaV, ((I_e_mean(i-1,j) + J_eex*params.Ex_Poisson_lambda)) * dt,...
                sqrt(I_e_var(i-1,j) + J_eex*params.Ex_Poisson_lambda) * sqrt(dt));
            
            %% 上一时刻的电流与电压卷积，计算每一个区间内神经元的跳跃
            e_newV = (V_lowlim+I_lowlim+dV):dV:(V_highlim+I_highlim);
            e_newV_p = conv(V_e_p, I_e_p);
            % 别忘了有可能会跳多个bin
            n_e_up(i,j) = round(n_e(i-1,j) * sum(e_newV_p(e_newV >= V_highlim))/sum(e_newV_p),digit_num);
            if n_e_up(i,j)
                V_e_up(i,j) = sum(e_newV(e_newV>=V_highlim) .* normal(e_newV_p(e_newV >= V_highlim)));
            end
            n_e_down(i,j) = round(n_e(i-1,j) * sum(e_newV_p(e_newV < V_lowlim))/sum(e_newV_p), digit_num);
            if n_e_down(i,j)
                V_e_down(i,j) = sum(e_newV(e_newV < V_lowlim) .* normal(e_newV_p(e_newV < V_lowlim)));
            end
            n_e_remain(i,j) = round(n_e(i-1,j) - n_e_down(i,j) - n_e_up(i,j),digit_num);     
            if n_e_remain(i,j)
                V_e_remain(i,j) = (n_e(i-1,j) * sum(e_newV .* normal(e_newV_p))...
                    - n_e_up(i,j) * V_e_up(i,j) ...
                    - n_e_down(i,j) * V_e_down(i,j))/n_e_remain(i,j);
                % 因为有四舍五入的环节，有可能因为这个导致留下的均值高于上限，需要额外处理
                if V_e_remain(i,j) >= V_highlim  
                    V_e_up(i,j) = (V_e_remain(i,j)*n_e_remain(i,j)+V_e_up(i,j)*n_e_up(i,j))...
                        /(n_e_remain(i,j)+n_e_up(i,j));
                    n_e_up(i,j) = n_e_remain(i,j)+n_e_up(i,j);
                    n_e_remain(i,j) = 0;
                elseif V_e_remain(i,j) < V_lowlim
                    V_e_down(i,j) = (V_e_remain(i,j)*n_e_remain(i,j) + V_e_down(i,j)*n_e_down(i,j))...
                        /(n_e_remain(i,j)+n_e_down(i,j));
                    n_e_down(i,j) = n_e_remain(i,j)+n_e_down(i,j);
                    n_e_remain(i,j) = 0;
                end
            end
        end
        
        %% I神经元
        if n_i(i-1,j) ~= 0 % 只有不为0的区间才有计算价值
            %% 计算V的概率分布（均匀分布）
            V_highlim = V_bin*(j+V_bin_min);
            V_lowlim = V_bin*(j+V_bin_min-1);
            if V_i_mean(i-1,j) < V_lowlim+dV
                V_i_p = zeros(1,(V_highlim-V_lowlim)/dV);
                V_i_p(1) = 1/dV;
            elseif V_i_mean(i-1,j) <= (V_highlim+V_lowlim)/2
                V_i_p = pdf('Uniform',(V_lowlim+dV):dV:V_highlim,V_lowlim,2*V_i_mean(i-1,j)-V_lowlim);
            elseif V_i_mean(i-1,j) >= V_highlim
                V_i_p = zeros(1,(V_highlim-V_lowlim)/dV);
                V_i_p(end) = 1/dV;
            else
                V_i_p = pdf('Uniform',(V_lowlim+dV):dV:V_highlim,2*V_i_mean(i-1,j)-V_highlim,V_highlim);
            end

            %% 计算上一时刻电流的分布参数（正态分布与泊松分布）
            dH_ie_mean = -H_ie_mean(i-1)/tau_ie + fre(i-1)*pie; % mean
            dH_ie_var = -H_ie_var(i-1)*2/(tau_ie) + H_ie_mean(i-1)/tau_ie + fre(i-1)*pie*(1-pie); % var
            H_ie_mean(i) = H_ie_mean(i-1) + dH_ie_mean*dt;
            H_ie_var(i) = H_ie_var(i-1) + dH_ie_var*dt;
            
            dH_ii_mean = -H_ii_mean(i-1)/tau_i + fri(i-1)*pii; % mean
            dH_ii_var = -H_ii_var(i-1)*2/(tau_i) + H_ii_mean(i-1)/tau_i + fri(i-1)*pii*(1-pii); % var
            H_ii_mean(i) = H_ii_mean(i-1) + dH_ii_mean*dt;
            H_ii_var(i) = H_ii_var(i-1) + dH_ii_var*dt;
            
            J_iex = 1;
            J_ie = params.s_ie;
            J_ii = (V_i_mean(i-1,j)+params.Mr)*params.s_ii/(params.M+params.Mr);    
            I_i_mean(i,j) = J_ie*H_ie_mean(i)/params.tau_ie - J_ii*H_ii_mean(i)/params.tau_i;
            I_i_var(i,j) = J_ie*H_ie_var(i)/params.tau_ie + J_ii*H_ii_var(i)/params.tau_i;
            I_i_p = pdf('Normal',IdeltaV,((I_i_mean(i-1,j) + J_iex*params.Ex_Poisson_lambda)) * dt,...
                sqrt(I_i_var(i-1,j) + J_iex*params.Ex_Poisson_lambda) * sqrt(dt));

            
            %% 计算每一个区间内神经元的跳跃
            i_newV = (V_lowlim+I_lowlim+dV):dV:(V_highlim+I_highlim);
            i_newV_p = conv(V_i_p, I_i_p);
            n_i_up(i,j) = round(n_i(i-1,j) * sum(i_newV_p(i_newV>=V_highlim))/sum(i_newV_p),digit_num);
            if n_i_up(i,j)
                V_i_up(i,j) = sum(i_newV(i_newV>=V_highlim) .* normal(i_newV_p(i_newV>=V_highlim)));
            end
            n_i_down(i,j) = round(n_i(i-1,j) * sum(i_newV_p(i_newV<V_lowlim))/sum(i_newV_p),digit_num);
            if n_i_down(i,j)
                V_i_down(i,j) = sum(i_newV(i_newV<V_lowlim) .* normal(i_newV_p(i_newV<V_lowlim)));
            end
            n_i_remain(i,j) = round(n_i(i-1,j) - n_i_down(i,j) - n_i_up(i,j),digit_num);     
            if n_i_remain(i,j)
                V_i_remain(i,j) = (n_i(i-1,j)*sum(i_newV .* normal(i_newV_p))...
                    - n_i_up(i,j)*V_i_up(i,j) ...
                    - n_i_down(i,j)*V_i_down(i,j))/n_i_remain(i,j);
                if V_i_remain(i,j)>=V_highlim
                    V_i_up(i,j) = (V_i_remain(i,j)*n_i_remain(i,j)+V_i_up(i,j)*n_i_up(i,j))/(n_i_remain(i,j)+n_i_up(i,j));
                    n_i_up(i,j) = n_i_remain(i,j)+n_i_up(i,j);
                    n_i_remain(i,j) = 0;
                elseif V_i_remain(i,j)<V_lowlim
                    V_i_down(i,j) = (V_i_remain(i,j)*n_i_remain(i,j)+V_i_down(i,j)*n_i_down(i,j))/(n_i_remain(i,j)+n_i_down(i,j));
                    n_i_down(i,j) = n_i_remain(i,j)+n_i_down(i,j);
                    n_i_remain(i,j) = 0;
                end
            end
        end
    end
    
    %% 计算跳跃后的结果
    for j = fliplr(1:(V_bin_num+1))
        if j == 1 % 第一个区间只有下降
            n_e(i,j) = n_e_remain(i,j) + n_e_down(i,j+1);
            V_e_mean(i,j) = (V_e_remain(i,j)*n_e_remain(i,j) + V_e_down(i,j+1)*n_e_down(i,j+1)) / n_e(i,j);
            n_i(i,j) = n_i_remain(i,j)  + n_i_down(i,j+1);
            V_i_mean(i,j) = (V_i_remain(i,j)*n_i_remain(i,j) + V_i_down(i,j+1)*n_i_down(i,j+1)) / n_i(i,j);
            
        elseif j == V_bin_num % 最后一个区间只有上升来的
            n_e(i,j) = n_e_remain(i,j) + n_e_up(i,j-1);
            V_e_mean(i,j) = (V_e_remain(i,j)*n_e_remain(i,j) + V_e_up(i,j-1)*n_e_up(i,j-1)) / n_e(i,j);
            n_i(i,j) = n_i_remain(i,j) + n_i_up(i,j-1);
            V_i_mean(i,j) = (V_i_remain(i,j)*n_i_remain(i,j) + V_i_up(i,j-1)*n_i_up(i,j-1)) / n_i(i,j);
            
        elseif j == V_bin_num+1 % 进入不应期的区间
            nfe(i) = n_e_up(i,j-1);
            fre(i) = nfe(i)/dt;
            ref_e(i) = ref_e(i-1) + nfe(i);
            nfi(i) = n_i_up(i,j-1);
            fri(i) = nfi(i)/dt;
            ref_i(i) = ref_i(i-1) + nfi(i); 
            
           %% 不应期神经元跳出
            if tau_r ~= 0
                dref_e = -ref_e(i)/tau_r + nfe(i);
                leave_e(i) = abs(round(dref_e*dt,digit_num));
                ref_e(i) = ref_e(i) - leave_e(i);

                dref_i = -ref_i(i)/tau_r + nfi(i);
                leave_i(i) = abs(round(dref_i*dt,digit_num));
                ref_i(i) = ref_i(i) - leave_i(i);
            else
                leave_e(i) = ref_e(i);
                ref_e(i) = 0;
                leave_i(i) = ref_i(i);
                ref_i(i) = 0;
            end
            
        elseif j == 1-V_bin_min % 0区间可能同时接受来自三个渠道的神经元
            n_e(i,j) = n_e_remain(i,j) + n_e_up(i,j-1) + n_e_down(i,j+1) + leave_e(i);
            V_e_mean(i,j) = (V_e_remain(i,j)*n_e_remain(i,j) ...
                + V_e_up(i,j-1)*n_e_up(i,j-1) ...
                + V_e_down(i,j+1)*n_e_down(i,j+1) ...
                + 0*leave_e(i))/n_e(i,j);
            n_i(i,j) = n_i_remain(i,j) + n_i_up(i,j-1) + n_i_down(i,j+1) + leave_i(i);
            V_i_mean(i,j) = (V_i_remain(i,j)*n_i_remain(i,j) ...
                + V_i_up(i,j-1)*n_i_up(i,j-1) ...
                + V_i_down(i,j+1)*n_i_down(i,j+1) ...
                + 0*leave_i(i))/n_i(i,j);

        else
            n_e(i,j) = n_e_remain(i,j) + n_e_up(i,j-1) + n_e_down(i,j+1);
            V_e_mean(i,j) = (V_e_remain(i,j)*n_e_remain(i,j) ...
                + V_e_up(i,j-1)*n_e_up(i,j-1) ...
                + V_e_down(i,j+1)*n_e_down(i,j+1))/n_e(i,j);
            n_i(i,j) = n_i_remain(i,j) + n_i_up(i,j-1) + n_i_down(i,j+1);
            V_i_mean(i,j) = (V_i_remain(i,j)*n_i_remain(i,j) ...
                + V_i_up(i,j-1)*n_i_up(i,j-1) ...
                + V_i_down(i,j+1)*n_i_down(i,j+1))/n_i(i,j);
        end
    end
    
end


res.n_e = n_e;
res.n_e_up = n_e_up;
res.n_e_down = n_e_down;
res.V_e_mean = V_e_mean;
res.V_e_all = V_e_mean.*n_e;
res.V_e_all(isnan(res.V_e_all)) = 0;
res.V_e_up = V_e_up;
res.V_e_remain = V_e_remain;
res.ref_e = ref_e;
res.H_ee_mean = H_ee_mean;
res.H_ei_mean = H_ei_mean;
res.H_ee_var = H_ee_var;
res.H_ei_var = H_ei_var;
res.nfe = nfe;
res.fre = fre;
J_eex = 1;
res.I_e_mean = I_e_mean+ J_iex*params.Ex_Poisson_lambda;
res.I_e_var = I_e_var+ J_iex*params.Ex_Poisson_lambda;

res.nfi = nfi;
res.fri = fri;
res.n_i = n_i;
res.V_i_mean = V_i_mean;
res.V_i_all = V_i_mean.*n_i;
res.V_i_all(isnan(res.V_i_all)) = 0;
res.ref_i = ref_i;
res.H_ie_mean = H_ie_mean;
res.H_ii_mean = H_ii_mean;
res.H_ie_var = H_ie_var;
res.H_ii_var = H_ii_var;
J_iex = 1;
res.I_i_mean = I_i_mean+ J_iex*params.Ex_Poisson_lambda;
res.I_i_var = I_i_var+ J_iex*params.Ex_Poisson_lambda;
end

function norm_list = normal(list)
norm_list = list/sum(list);
end
