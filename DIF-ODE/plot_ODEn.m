x = 5;
if x == 1
    V_bin = 2;
    V_bin_min = -3;
    V_bin_num = params.M/V_bin - V_bin_min;
    t = size(res_lif.VE);
    n_e = zeros(t(1),V_bin_num+1);
    n_i = zeros(t(1),V_bin_num+1);
    for i = 1:t(1)
       n_e(i,1:V_bin_num) = histcounts(res_lif.VE(i,:),V_bin_min*V_bin:V_bin:params.M); 
       n_e(i,V_bin_num+1) = sum(res_lif.VE(i,:)==0);
       n_e(i,1-V_bin_min) = n_e(i,1-V_bin_min) - n_e(i,V_bin_num+1);
       n_i(i,1:V_bin_num) = histcounts(res_lif.VI(i,:),V_bin_min*V_bin:V_bin:params.M); 
       n_i(i,V_bin_num+1) = sum(res_lif.VI(i,:)==0);
       n_i(i,1-V_bin_min) = n_i(i,1-V_bin_min) - n_i(i,V_bin_num+1);
    end

elseif x==2
    V_bin = 10;
    V_bin_min = -3;
    V_bin_num = params.M/V_bin - V_bin_min;
    t = size(res_Ifull_model.V_e);
    n_e = zeros(t(1),V_bin_num+1);
    n_i = zeros(t(1),V_bin_num+1);
    for i = 1:t(1)
       n_e(i,1:V_bin_num) = histcounts(res_Ifull_model.V_e(i,:),V_bin_min*V_bin:V_bin:params.M); 
       n_e(i,V_bin_num+1) = sum(res_Ifull_model.V_e(i,:)==0);
       n_e(i,1-V_bin_min) = n_e(i,1-V_bin_min) - n_e(i,V_bin_num+1);
       n_i(i,1:V_bin_num) = histcounts(res_Ifull_model.V_i(i,:),V_bin_min*V_bin:V_bin:params.M); 
       n_i(i,V_bin_num+1) = sum(res_Ifull_model.V_i(i,:)==0);
       n_i(i,1-V_bin_min) = n_i(i,1-V_bin_min) - n_i(i,V_bin_num+1);
    end
elseif x == 3
    V_bin = 10;
    V_bin_min = -4;
    V_bin_num = params.M/V_bin - V_bin_min;
    t = size(res_Ifull_model.V_e);
    n_e = zeros(t(1),V_bin_num+1);
    n_i = zeros(t(1),V_bin_num+1);
    for i = 1:t(1)
       n_e(i,1:V_bin_num) = histcounts(res_Ifull_model.V_e(i,:),V_bin_min*V_bin:V_bin:params.M); 
       n_e(i,V_bin_num+1) = sum(res_Ifull_model.V_e(i,:)==0);
       n_e(i,1-V_bin_min) = n_e(i,1-V_bin_min) - n_e(i,V_bin_num+1);
       n_i(i,1:V_bin_num) = histcounts(res_Ifull_model.V_i(i,:),V_bin_min*V_bin:V_bin:params.M); 
       n_i(i,V_bin_num+1) = sum(res_Ifull_model.V_i(i,:)==0);
       n_i(i,1-V_bin_min) = n_i(i,1-V_bin_min) - n_i(i,V_bin_num+1);
    end
    n_e2 = [res_ODEfull_modelhom.n_e res_ODEfull_modelhom.ref_e];
    n_i2 = [res_ODEfull_modelhom.n_i res_ODEfull_modelhom.ref_i];
elseif x==4
    n_e = [res_ODEfull_model.n_e res_ODEfull_model.ref_e];
    n_i = [res_ODEfull_model.n_i res_ODEfull_model.ref_i];
    n_e2 = [res_ODEfull_model2.n_e res_ODEfull_model2.ref_e];
    n_i2 = [res_ODEfull_model2.n_i res_ODEfull_model2.ref_i];
else
    n_e = [res_DIFODE.n_e res_DIFODE.ref_e];
    n_i = [res_DIFODE.n_i res_DIFODE.ref_i];
end


% n_e = [res_ODEfull_modelhomroundnofixdt01.n_e res_ODEfull_modelhomroundnofixdt01.ref_e];
% n_i = [res_ODEfull_modelhomroundnofixdt01.n_i res_ODEfull_modelhomroundnofixdt01.ref_i];
% X = categorical({'[-20,0)','[0,20)','[20,40)','[40,60)','[60,80)','[80,100)','ref'});
% X = reordercats(X,{'[-20,0)','[0,20)','[20,40)','[40,60)','[60,80)','[80,100)','ref'});
X = categorical({'[-40,-30)','[-30,-20)','[-20,-10)','[-10,0)','[0,10)','[10,20)','[20,30)','[30,40)','[40,50)','[50,60)','[60,70)','[70,80)','[80,90)','[90,100)','ref'});
X = reordercats(X,{'[-40,-30)','[-30,-20)','[-20,-10)','[-10,0)','[0,10)','[10,20)','[20,30)','[30,40)','[40,50)','[50,60)','[60,70)','[70,80)','[80,90)','[90,100)','ref'});
% X = categorical({'[-12,-8)','[-8,-4)','[-4,0)','[0,4)','[4,8)','[8,12)','[12,16)','[16,20)',...
%     '[20,24)','[24,28)','[28,32)','[32,36)','[36,40)','[40,44)','[44,48)','[48,52)','[52,56)','[56,60)',...
%     '[60,64)','[64,68)','[68,72)','[72,76)','[76,80)','[80,84)','[84,88)','[88,92)','[92,96)','[96,100)','ref'});
% X = reordercats(X,{'[-12,-8)','[-8,-4)','[-4,0)','[0,4)','[4,8)','[8,12)','[12,16)','[16,20)',...
%     '[20,24)','[24,28)','[28,32)','[32,36)','[36,40)','[40,44)','[44,48)','[48,52)','[52,56)','[56,60)',...
%     '[60,64)','[64,68)','[68,72)','[72,76)','[76,80)','[80,84)','[84,88)','[88,92)','[92,96)','[96,100)','ref'});
figure
xlabel('V bin')
ylabel('neuron number')
dt = params.dt;
% dt2 = 0.08;
step = 10;
start_time = 500;
for i = (start_time/dt):step:((start_time+2000)/dt)
    Y1 = n_e(i,:);
    subplot(2,2,1)
    b = bar(X,Y1);
    b.FaceColor = 'flat';
    b.CData(end,:) = [.5 0 .2];
    ylim([0 params.ne])
    title('ODE N E')
    subplot(2,2,2)
    Y2 = n_i(i,:);
    b = bar(X,Y2);
    b.FaceColor = 'flat';
    b.CData(end,:) = [.5 0 .2];
    title('ODE N I')
    ylim([0 params.ni])
    
    
    
%     Y1 = n_e2(i,:);
%     subplot(2,2,3)
%     b = bar(X,Y1);
%     b.FaceColor = 'flat';
%     b.CData(end,:) = [.5 0 .2];
%     ylim([0 300])
%     title('ODE N E')
%     subplot(2,2,4)
%     Y2 = n_i2(i,:);
%     b = bar(X,Y2);
%     b.FaceColor = 'flat';
%     b.CData(end,:) = [.5 0 .2];
%     title('ODE N I')
%     ylim([0 100])
%     
    
    sgtitle(strcat('time',num2str(i*dt)));


pause(0.01);
end
