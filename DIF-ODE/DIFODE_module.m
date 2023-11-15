function [output] = DIFODE_module(state,params)
output = struct;

%% preparation
dt = params.dt;
tau_r = params.tau_r;
V_bin = params.V_bin;
V_bin_min = params.V_bin_min;
V_bin_num = params.V_bin_num;
digit_num = params.digit_num;
t_end = 2;
i = 2;
V_n_mean = zeros(t_end, V_bin_num);

V_n_remain = zeros(t_end, V_bin_num);
V_n_up = zeros(t_end, V_bin_num);
V_n_down = zeros(t_end, V_bin_num);
I_n_mean = zeros(t_end, V_bin_num);
I_n_var = zeros(t_end, V_bin_num);


n_n = zeros(t_end, V_bin_num);
n_n_remain = zeros(t_end, V_bin_num);
n_n_up = zeros(t_end, V_bin_num);
n_n_down = zeros(t_end, V_bin_num);

nf_n = zeros(t_end, 1);
fr_n = zeros(t_end, 1);
ref_n = zeros(t_end, 1);
leave_n = zeros(t_end, 1);


n_n(1,:) = state.n_n;
V_n_all = state.V_n_all;
V_n_mean(1,:) = V_n_all./n_n(1,:);
fr_n(1) = state.fr_n(1);
nf_n(1) = state.fr_n(1)*params.dt;
ref_n(1) = state.ref_n;
I_n_mean(1,:) = state.I_n_mean; 
I_n_var(1,:) = state.I_n_var;

%% run model
%% 计算每个区间神经元电流
for j = 1:V_bin_num
    if n_n(i-1,j) ~= 0 % 只有不为0的区间才有计算价值
        %% 计算V的概率分布（均匀分布）
        V_highlim = V_bin*(j+V_bin_min);
        V_lowlim = V_bin*(j+V_bin_min-1);
        if V_n_mean(i-1,j) < V_lowlim+0.0001
            a = V_lowlim;
            b = V_lowlim+0.0001;
        elseif V_n_mean(i-1,j) <= (V_highlim+V_lowlim)/2
            a = V_lowlim;
            b = 2*V_n_mean(i-1,j)-V_lowlim;
        elseif V_n_mean(i-1,j) >= V_highlim
            a = V_highlim-0.0001;
            b = V_highlim;
        else
            a = 2*V_n_mean(i-1,j)-V_highlim;
            b = V_highlim;
        end
        mu = I_n_mean(1,j) * dt;
        sigma = sqrt(I_n_var(1,j)) * sqrt(dt);
        %% 上一时刻的电流与电压卷积，计算每一个区间内神经元的跳跃
        % 别忘了有可能会跳多个bin，计算完成后要检查一下max(res_DIFODE.I_e_mean*dt)
        up_probability_accumulation = (0.5-(sqrt(2/pi)*sigma*(exp(-(b+mu-V_highlim)^2./(2*sigma^2))-exp(-(a+mu-V_highlim)^2./(2.*sigma^2)))-...
            (a+mu-V_highlim)*erf((a+mu-V_highlim)/(sqrt(2)*sigma))+(b+mu-V_highlim)*erf((b+mu-V_highlim)/(sqrt(2)*sigma)))/(2*(a-b)));
        if abs(up_probability_accumulation) < 10^(-15)% 因为数值计算的关系，有时候0.5-0.5会出现一个极小的数字，要去掉
            up_probability_accumulation = 0;
        end
        n_n_up(i,j) = round(n_n(i-1,j)*up_probability_accumulation,digit_num);
        if n_n_up(i,j)
            upV_probability_accumulation = (mu+(a+b)/2)/2-(a^2*erf((-a-mu+V_highlim)/(sqrt(2)*sigma))+mu^2*erf((-a-mu+V_highlim)/(sqrt(2)*sigma))...
                +sigma^2*erf((-a-mu+V_highlim)/(sqrt(2)*sigma))-sqrt(2/pi)*sigma*(a+mu+V_highlim)*exp(-(a+mu-V_highlim)^2/(2*sigma^2))...
                +V_highlim^2*erf((a+mu-V_highlim)/(sqrt(2)*sigma))+2*a*mu*erf((-a-mu+V_highlim)/(sqrt(2)*sigma))...
                -b^2.*erf((-b-mu+V_highlim)/(sqrt(2)*sigma))-mu^2*erf((-b-mu+V_highlim)/(sqrt(2)*sigma))...
                -sigma^2*erf((-b-mu+V_highlim)/(sqrt(2)*sigma))+sqrt(2/pi)*sigma*(b+mu+V_highlim)*exp(-(b+mu-V_highlim)^2/(2*sigma^2))...
                -V_highlim^2*erf((b+mu-V_highlim)/(sqrt(2)*sigma))-2*b*mu*erf((-b-mu+V_highlim)/(sqrt(2)*sigma)))/(4*(a-b));
            V_n_up(i,j) = upV_probability_accumulation/up_probability_accumulation;
        else
            upV_probability_accumulation = 0;
        end
        down_probability_accumulation=((sqrt(2/pi)*sigma*(exp(-(b+mu-V_lowlim)^2./(2*sigma^2))-exp(-(a+mu-V_lowlim)^2./(2.*sigma^2)))-...
            (a+mu-V_lowlim)*erf((a+mu-V_lowlim)/(sqrt(2)*sigma))+(b+mu-V_lowlim)*erf((b+mu-V_lowlim)/(sqrt(2)*sigma)))/(2*(a-b))+0.5);
        if abs(down_probability_accumulation) < 10^(-15)% 因为数值计算的关系，有时候0.5-0.5会出现一个极小的数字，要去掉
            down_probability_accumulation = 0;
        end
        n_n_down(i,j) = round(n_n(i-1,j)*down_probability_accumulation,digit_num);
        if n_n_down(i,j)
            downV_probability_accumulation = (a^2*erf((-a-mu+V_lowlim)/(sqrt(2)*sigma))+mu^2*erf((-a-mu+V_lowlim)/(sqrt(2)*sigma))...
                +sigma^2*erf((-a-mu+V_lowlim)/(sqrt(2)*sigma))-sqrt(2/pi)*sigma*(a+mu+V_lowlim)*exp(-(a+mu-V_lowlim)^2/(2*sigma^2))...
                +V_lowlim^2*erf((a+mu-V_lowlim)/(sqrt(2)*sigma))+2*a*mu*erf((-a-mu+V_lowlim)/(sqrt(2)*sigma))...
                -b^2.*erf((-b-mu+V_lowlim)/(sqrt(2)*sigma))-mu^2*erf((-b-mu+V_lowlim)/(sqrt(2)*sigma))...
                -sigma^2*erf((-b-mu+V_lowlim)/(sqrt(2)*sigma))+sqrt(2/pi)*sigma*(b+mu+V_lowlim)*exp(-(b+mu-V_lowlim)^2/(2*sigma^2))...
                -V_lowlim^2*erf((b+mu-V_lowlim)/(sqrt(2)*sigma))-2*b*mu*erf((-b-mu+V_lowlim)/(sqrt(2)*sigma)))/(4*(a-b))+(mu+(a+b)/2)/2;
            V_n_down(i,j) = downV_probability_accumulation/down_probability_accumulation;
        else
            downV_probability_accumulation = 0;
        end
        n_n_remain(i,j) = round(n_n(i-1,j) - n_n_down(i,j) - n_n_up(i,j),digit_num);
        remain_probability_accumulation = n_n_remain(i,j)/n_n(i-1,j);
        if n_n_remain(i,j)
            remainV_probability_accumulation = (mu+(a+b)/2)-upV_probability_accumulation-downV_probability_accumulation;
            V_n_remain(i,j) = remainV_probability_accumulation/remain_probability_accumulation;
            %% 因为四舍五入的问题，有时候会导致留下的电压大于最高值，这一步是为了把多余的电压分摊过去
            if V_n_remain(i,j) >= V_highlim
                V_n_up(i,j) = (V_n_remain(i,j)*n_n_remain(i,j)+V_n_up(i,j)*n_n_up(i,j))...
                    /(n_n_remain(i,j)+n_n_up(i,j));
                n_n_up(i,j) = n_n_remain(i,j)+n_n_up(i,j);
                n_n_remain(i,j) = 0;
            elseif V_n_remain(i,j) < V_lowlim
                V_n_down(i,j) = (V_n_remain(i,j)*n_n_remain(i,j) + V_n_down(i,j)*n_n_down(i,j))...
                    /(n_n_remain(i,j)+n_n_down(i,j));
                n_n_down(i,j) = n_n_remain(i,j)+n_n_down(i,j);
                n_n_remain(i,j) = 0;
            end
        end
    end
end


%% 计算跳跃后的结果
for j = fliplr(1:(V_bin_num+1))
    if j == 1 % 第一个区间只有下降
        n_n(i,j) = n_n_remain(i,j) + n_n_down(i,j+1);
        V_n_mean(i,j) = (V_n_remain(i,j)*n_n_remain(i,j) + V_n_down(i,j+1)*n_n_down(i,j+1)) / n_n(i,j);

    elseif j == V_bin_num % 最后一个区间只有上升来的
        n_n(i,j) = n_n_remain(i,j) + n_n_up(i,j-1);
        V_n_mean(i,j) = (V_n_remain(i,j)*n_n_remain(i,j) + V_n_up(i,j-1)*n_n_up(i,j-1)) / n_n(i,j);

    elseif j == V_bin_num+1 % 进入不应期的区间
        nf_n(i) = n_n_up(i,j-1);
        fr_n(i) = nf_n(i)/dt;
        ref_n(i) = ref_n(i-1) + nf_n(i);

       %% 不应期神经元跳出
        if tau_r ~= 0
            dref = -ref_n(i)/tau_r + nf_n(i);
            leave_n(i) = abs(round(dref*dt,digit_num));
            ref_n(i) = ref_n(i) - leave_n(i);
        else
            leave_n(i) = ref_n(i);
            ref_n(i) = 0;
        end

    elseif j == 1-V_bin_min % 0区间可能同时接受来自三个渠道的神经元
        n_n(i,j) = n_n_remain(i,j) + n_n_up(i,j-1) + n_n_down(i,j+1) + leave_n(i);
        V_n_mean(i,j) = (V_n_remain(i,j)*n_n_remain(i,j) ...
            + V_n_up(i,j-1)*n_n_up(i,j-1) ...
            + V_n_down(i,j+1)*n_n_down(i,j+1) ...
            + 0*leave_n(i))/n_n(i,j);
    else
        n_n(i,j) = n_n_remain(i,j) + n_n_up(i,j-1) + n_n_down(i,j+1);
        V_n_mean(i,j) = (V_n_remain(i,j)*n_n_remain(i,j) ...
            + V_n_up(i,j-1)*n_n_up(i,j-1) ...
            + V_n_down(i,j+1)*n_n_down(i,j+1))/n_n(i,j);
    end
end


output.n_n = n_n(2,:);
output.V_n_mean = V_n_mean(2,:);
output.V_n_all = V_n_mean(2,:).*n_n(2,:);
output.V_n_all(isnan(output.V_n_all)) = 0;
output.ref_n = ref_n(2);
output.nf_n = nf_n(2);
output.fr_n = fr_n(2);
output.I_n_mean = I_n_mean(2);
output.I_n_var = I_n_var(2);
end
