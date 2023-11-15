x = 4;
start = struct;
time = 400;
dt = 0.01;
i = time/dt;
V_bin = 10;
V_bin_min = -4;
V_bin_num = params.M/V_bin - V_bin_min;
start.n_e = zeros(1,V_bin_num);
start.n_i = zeros(1,V_bin_num);
start.V_e_mean = zeros(1,V_bin_num);
start.V_i_mean = zeros(1,V_bin_num);

if x == 1 % 1代表LIF1
    fr = firing_rate(res_lif,param3,'part');
    start.nfe = fr.e(i);
    start.nfi = fr.i(i);
    start.ref_e = 0;
    start.ref_i = 0;
    for j = 1:V_bin_num
        V_highlim = V_bin*(j+V_bin_min);
        V_lowlim = V_bin*(j+V_bin_min-1);
        Ejbin_index =  res_lif.VE(i,:) < V_highlim & res_lif.VE(i,:) >= V_lowlim;
        start.n_e(j) = sum(Ejbin_index);
        start.V_e_mean(j) = mean(res_lif.VE(i,Ejbin_index));
        Ijbin_index =  res_lif.VI(i,:) < V_highlim & res_lif.VI(i,:) >= V_lowlim;
        start.n_i(j) = sum(Ijbin_index);
        start.V_i_mean(j) = mean(res_lif.VI(i,Ijbin_index));
    end

    start.n_e(1-V_bin_min) = start.n_e(1-V_bin_min) - start.ref_e;
    start.n_i(1-V_bin_min) = start.n_i(1-V_bin_min) - start.ref_i;
    start.V_e_all = start.n_e.*start.V_e_mean;
    start.V_i_all = start.n_i.*start.V_i_mean;
    start.H_ee_mean = mean(res_lif.H_ee(i,:));
    start.H_ei_mean = mean(res_lif.H_ei(i,:));
    start.H_ie_mean = mean(res_lif.H_ie(i,:));
    start.H_ii_mean = mean(res_lif.H_ii(i,:));
    start.H_ee_var = var(res_lif.H_ee(i,:),0,2);
    start.H_ei_var = var(res_lif.H_ei(i,:),0,2);
    start.H_ie_var = var(res_lif.H_ie(i,:),0,2);
    start.H_ii_var = var(res_lif.H_ii(i,:),0,2);

    
elseif x == 2
    start.nfe = res_Ifull_model.nfe(i);
    start.nfi = res_Ifull_model.nfi(i);
    start.ref_e = res_Ifull_model.ref_e(i);
    start.ref_i = res_Ifull_model.ref_i(i);
    for j = 1:V_bin_num

        V_highlim = V_bin*(j+V_bin_min);
        V_lowlim = V_bin*(j+V_bin_min-1);
        Ejbin_index =  res_Ifull_model.V_e(i,:) < V_highlim & res_Ifull_model.V_e(i,:) >= V_lowlim;
        start.n_e(j) = sum(Ejbin_index);
        start.V_e_mean(j) = mean(res_Ifull_model.V_e(i,Ejbin_index));


        Ijbin_index =  res_Ifull_model.V_i(i,:) < V_highlim & res_Ifull_model.V_i(i,:) >= V_lowlim;
        start.n_i(j) = sum(Ijbin_index);
        start.V_i_mean(j) = mean(res_Ifull_model.V_i(i,Ijbin_index));

    end
    start.n_e(1-V_bin_min) = start.n_e(1-V_bin_min) - start.ref_e;
    start.n_i(1-V_bin_min) = start.n_i(1-V_bin_min) - start.ref_i;
    start.V_e = res_Ifull_model.V_e(i,:);
    start.V_i = res_Ifull_model.V_i(i,:);
    start.H_ee = mean(res_Ifull_model.H_ee(i,:));
    start.H_ei = mean(res_Ifull_model.H_ei(i,:));
    start.H_ie = mean(res_Ifull_model.H_ie(i,:));
    start.H_ii = mean(res_Ifull_model.H_ii(i,:));
    start.H_ee_var = var(res_Ifull_model.H_ee(i,:),0,2);
    start.H_ei_var = var(res_Ifull_model.H_ei(i,:),0,2);
    start.H_ie_var = var(res_Ifull_model.H_ie(i,:),0,2);
    start.H_ii_var = var(res_Ifull_model.H_ii(i,:),0,2);
    % start.ref_e = sum(res_Ifull_model.V_e(i,:) == 0);
    % start.ref_i = sum(res_Ifull_model.V_i(i,:) == 0);
    
    
elseif x == 3
    start.H_ee_mean = res_ODEfull_model.H_ee_mean(i);
    start.H_ei_mean = res_ODEfull_model.H_ei_mean(i);
    start.H_ie_mean = res_ODEfull_model.H_ie_mean(i);
    start.H_ii_mean = res_ODEfull_model.H_ii_mean(i);
    start.H_ee_var = res_ODEfull_model.H_ee_var(i);
    start.H_ei_var = res_ODEfull_model.H_ei_var(i);
    start.H_ie_var = res_ODEfull_model.H_ie_var(i);
    start.H_ii_var = res_ODEfull_model.H_ii_var(i);
    
    start.n_e = res_ODEfull_model.n_e(i,:);
    start.n_i = res_ODEfull_model.n_i(i,:);
    start.V_e_all = res_ODEfull_model.V_e_all(i,:);
    start.V_i_all = res_ODEfull_model.V_i_all(i,:);
    start.V_e_mean = res_ODEfull_model.V_e_mean(i,:);
    start.V_i_mean = res_ODEfull_model.V_i_mean(i,:);
    start.fre = res_ODEfull_model.fre(i);
    start.fri = res_ODEfull_model.fri(i);
    start.nfe = res_ODEfull_model.nfe(i);
    start.nfi = res_ODEfull_model.nfi(i);
    start.ref_e = res_ODEfull_model.ref_e(i);
    start.ref_i = res_ODEfull_model.ref_i(i);
else
    start.H_ee_mean = res_DIFODE.H_ee_mean(i);
    start.H_ei_mean = res_DIFODE.H_ei_mean(i);
    start.H_ie_mean = res_DIFODE.H_ie_mean(i);
    start.H_ii_mean = res_DIFODE.H_ii_mean(i);
    start.H_ee_var = res_DIFODE.H_ee_var(i);
    start.H_ei_var = res_DIFODE.H_ei_var(i);
    start.H_ie_var = res_DIFODE.H_ie_var(i);
    start.H_ii_var = res_DIFODE.H_ii_var(i);
    
    start.n_e = res_DIFODE.n_e(i,:);
    start.n_i = res_DIFODE.n_i(i,:);
    start.V_e_all = res_DIFODE.V_e_all(i,:);
    start.V_i_all = res_DIFODE.V_i_all(i,:);
    start.V_e_mean = res_DIFODE.V_e_mean(i,:);
    start.V_i_mean = res_DIFODE.V_i_mean(i,:);
    start.fr_e = res_DIFODE.fr_e(i);
    start.fr_i = res_DIFODE.fr_i(i);
    start.ref_e = res_DIFODE.ref_e(i);
    start.ref_i = res_DIFODE.ref_i(i);
end



