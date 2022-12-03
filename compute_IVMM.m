function [estate_mat,istate_mat] = compute_IVMM(res,params,estate_mat_before,istate_mat_before)
if nargin > 4
    error('输入变量过多！');
elseif nargin == 2
    flag = 0; % 默认情况下无before
else 
    flag = 1;
end

%% 参数
ne = params.ne;
ni = params.ni;
tindex_end = params.duration_time/0.1;
all_time = res.t(1:tindex_end);
tau_ee =  params.tau_ee;
tau_ie = params.tau_ie;
tau_i = params.tau_i;
switch res.Ex_type
    case{'constant'}
        tau_ex = 1/params.Ex_Poisson_lambda;
    case{'sin'}
        C1 = params.Ex_sin_C1;
        C2 = params.Ex_sin_C2;
        C3 = params.Ex_sin_C3;
        C4 = params.Ex_sin_C4;
        lambda = C1*sin(C2*all_time+C3)+C4; % Time of all maximum points
        tau_ex = 1./lambda';
end


%% 计算ref
fr = firing_rate(res,params,'part');
nfe = fr.e;
nfi = fr.i;
if isfield(res, 'iref')
    eref = reshape(res.eref, ne, length(res.eref)/ne);
    eref = eref';
    iref = reshape(res.iref, ni, length(res.iref)/ni);
    iref = iref';
    ref_e = mean(eref,2,'omitnan');
    ref_i = mean(iref,2,'omitnan');
else
    tau_ref = 3;
    ht = 0.1:0.1:10;
    ht2 = 0:0.1:9.9;
    h_tau_ref = (exp(-(1/tau_ref)*ht)+exp(-(1/tau_ref)*ht2))/2;
    ref_e = conv(nfe,h_tau_ref);
    ref_e = ref_e(1:end-length(ht)+1)/ne;
    ref_i = conv(nfi,h_tau_ref);
    ref_i = ref_i(1:end-length(ht)+1)/ni;
end


%% H current
H_ee = res.H_ee(1:tindex_end*ne);
H_ie = res.H_ie(1:tindex_end*ni);
H_ei = res.H_ei(1:tindex_end*ne);
H_ii = res.H_ii(1:tindex_end*ni);
V_e = res.V_e(1:tindex_end*ne);
V_i = res.V_i(1:tindex_end*ni);

J_eex = ones(1,length(all_time));
J_iex = ones(1,length(all_time));
J_ee = params.s_ee;
J_ie = params.s_ie;
J_ei = (V_e+params.Mr)*params.s_ei/(params.M+params.Mr);
J_ii = (V_i+params.Mr)*params.s_ii/(params.M+params.Mr);

I_eex = J_eex'./tau_ex;
I_iex = J_iex'./tau_ex;
I_ee = J_ee.*H_ee/tau_ee;
I_ie = J_ie.*H_ie/tau_ie;
I_ei = J_ei.*H_ei/tau_i;
I_ii = J_ii.*H_ii/tau_i;

V_e = reshape(V_e, ne, length(V_e)/ne);
V_e = V_e';
V_i = reshape(V_i, ni, length(V_i)/ni);
V_i = V_i';
I_ee = reshape(I_ee, ne, length(I_ee)/ne);
I_ee = I_ee';
I_ie = reshape(I_ie, ni, length(I_ie)/ni);
I_ie = I_ie';
I_ei = reshape(I_ei, ne, length(I_ei)/ne);
I_ei = I_ei';
I_ii = reshape(I_ii, ni, length(I_ii)/ni);
I_ii = I_ii';
V_e(V_e==0) = NaN;
V_i(V_i==0) = NaN;

V_e = mean(V_e,2,'omitnan');
V_i = mean(V_i,2,'omitnan');
I_ee = mean(I_ee,2);
I_ie = mean(I_ie,2);
I_ei = mean(I_ei,2);
I_ii = mean(I_ii,2);

I_e = I_ee+I_eex-I_ei;
I_i = I_ie+I_iex-I_ii;

%% real current
% V_e = res.V_e(1:tindex_end*ne);
% V_i = res.V_i(1:tindex_end*ni);
% V_e = reshape(V_e, ne, length(V_e)/ne);
% V_e = V_e';
% V_i = reshape(V_i, ni, length(V_i)/ni);
% V_i = V_i';
% V_e_var = std(V_e,0,2,'omitnan');
% V_e = mean(V_e,2,'omitnan');
% V_i_var = std(V_i,0,2,'omitnan');
% V_i = mean(V_i,2,'omitnan');
% I_eex = 1*histcounts(res.t1_Exsti(2:end,1:ne),[0,all_time(1:tindex_end)])./diff([0,all_time(1:tindex_end)]);
% I_iex = 1*histcounts(res.t1_Exsti(2:end,ne+1:ne+ni),[0,all_time(1:tindex_end)])./diff([0,all_time(1:tindex_end)]);
% I_ee = params.s_ee*histcounts(res.t2_Esti(2:end,1:ne),[0,all_time(1:tindex_end)])./diff([0,all_time(1:tindex_end)]);
% I_ie = params.s_ie*histcounts(res.t2_Esti(2:end,ne+1:ne+ni),[0,all_time(1:tindex_end)])./diff([0,all_time(1:tindex_end)]); 
% J_ei = (V_e+params.Mr)*params.s_ei/(params.M+params.Mr);
% J_ii = (V_i+params.Mr)*params.s_ii/(params.M+params.Mr);
% I_ei = J_ei'.*histcounts(res.t3_Isti(2:end,1:ne),[0,all_time(1:tindex_end)])./diff([0,all_time(1:tindex_end)]);
% I_ii = J_ii'.*histcounts(res.t3_Isti(2:end,ne+1:ne+ni),[0,all_time(1:tindex_end)])./diff([0,all_time(1:tindex_end)]); 
% I_e = (I_ee+I_eex-I_ei)/ne;
% I_i = (I_ie+I_iex-I_ii)/ni;


% I_e_noref = I_e(:) .* (1-ref_e(:));
% I_i_noref = I_i(:) .* (1-ref_i(:));
I_e_noref = I_e(:);
I_i_noref = I_e(:);

%% 计算转移情况（分了两种，一种同时统计V和nf，一种只统计V
V_e = round(V_e);
V_i = round(V_i);
I_e_noref = round(I_e_noref);
I_i_noref = round(I_i_noref);
min_V_e = -20;
max_V_e = 100;
min_V_i = -20;
max_V_i = 100;
min_I_e_noref = -30;
min_I_i_noref = -30;
max_I_e_noref = 120;
max_I_i_noref = 120;
V_e_statenum = max_V_e-min_V_e+1;
V_i_statenum = max_V_i-min_V_i+1;
I_e_noref_statenum = max_I_e_noref-min_I_e_noref+1;
I_i_noref_statenum = max_I_i_noref-min_I_i_noref+1;

% 统计E态的数目
if flag
    estate_mat = estate_mat_before;
else
    estate_mat = zeros(V_e_statenum,I_e_noref_statenum,V_e_statenum);
end
% 只统计Vi+1
for i = 300:length(V_e)-1
    estate_mat(V_e(i)-min_V_e+1,I_e_noref(i)-min_I_e_noref+1,V_e(i+1)-min_V_e+1) = ...
        estate_mat(V_e(i)-min_V_e+1,I_e_noref(i)-min_I_e_noref+1,V_e(i+1)-min_V_e+1)+1;
end

% 统计I态的数目
if flag
    istate_mat = istate_mat_before;
else
    istate_mat = zeros(V_i_statenum,I_i_noref_statenum,V_i_statenum);
end

% 只统计Vi+1
for i = 300:length(V_e)-1
    istate_mat(V_i(i)-min_V_i+1,I_i_noref(i)-min_I_i_noref+1,V_i(i+1)-min_V_i+1) = ...
        istate_mat(V_i(i)-min_V_i+1,I_i_noref(i)-min_I_i_noref+1,V_i(i+1)-min_V_i+1)+1;
end

%% 合并同类项（只统计V就不需要这一步了）
% if flag
%     estate_mat = estate_mat_before;
% else
%     estate_mat = cell(V_e_statenum,I_e_noref_statenum);
% end
% if flag
%     istate_mat = istate_mat_before;
% else
%     istate_mat = cell(V_i_statenum,I_i_noref_statenum);
% end
% % 同时统计Vi+1与Nfi+1
% for i = 300:length(V_e)-1
%     estate_mat{V_e(i)-min_V_e+1,I_e_noref(i)-min_I_e_noref+1} = [estate_mat{V_e(i)-min_V_e+1,I_e_noref(i)-min_I_e_noref+1};[V_e(i+1) nfe(i+1) 1]];
% end
% % 同时统计Vi+1与Nfi+1
% for i = 300:length(V_e)-1
%     istate_mat{V_i(i)-min_V_i+1,I_i_noref(i)-min_I_i_noref+1} = [istate_mat{V_i(i)-min_V_i+1,I_i_noref(i)-min_I_i_noref+1};[V_i(i+1) nfi(i+1) 1]];
% end

% for i = 1:V_e_statenum*I_e_noref_statenum
%     estate = estate_mat{i};
%     if ~isempty(estate)
%         [unique_rows,~,ind] = unique(estate(:,1:2),'rows');
%         ecounts = [];
%         for j = 1:length(unique(ind))
%             ecounts = [ecounts; sum(estate(ind==j,3))];
%         end
%         arr = [unique_rows,ecounts];
%         [temp,~]=sortrows(arr,[1 -3]);
%         estate_mat{i} = temp;
%     end
% end

% for i = 1:V_i_statenum*I_i_noref_statenum
%     istate = istate_mat{i};
%     if ~isempty(istate)
%         [unique_rows,~,ind] = unique(istate(:,1:2),'rows');
%         icounts = [];
%         for j = 1:length(unique(ind))
%             icounts = [icounts; sum(istate(ind==j,3))];
%         end
%         arr = [unique_rows,icounts];
%         [temp,~]=sortrows(arr,[1 -3]);
%         istate_mat{i} = temp;
%     end
% end

end
