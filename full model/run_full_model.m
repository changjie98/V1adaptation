function [res] = run_full_model(params, need_set, res_before)
if nargin > 3
    error('输入变量过多！');
elseif nargin == 3
    havebefore = true;
else
    havebefore = false;
end

%% preparation

t = 0;
res.V_e = [];
res.V_i = [];
res.H_ie = [];
res.H_ii = [];
res.H_ee = [];
res.H_ei = [];
res.t = [];
res.eref = [];
res.iref = [];


c = zeros(4,params.ne + params.ni);
% c(1,:) = Ex_input_tau(t, "sin");
% c(1,:) = repmat(1/params.Ex_Poisson_lambda, 1, params.ne + params.ni);

c(2,1:params.ne)                   = params.tau_ee;
c(2,params.ne+1:params.ne+params.ni) = params.tau_ie;
c(3,:)                            = params.tau_i;
c(4,:)                            = params.tau_r;
duration_time                     = params.duration_time;
%c is time constant matrix. The first row is lambda_e(i), the second row
%correspond to H_e, the third row correspond to H_i, the fourth row
%correspond to tau_r.

if havebefore
    s = zeros(4,params.ne+params.ni);
    s(1,:) = [res_before.V_e(end-params.ne+1,end) res_before.V_i(end-params.ni+1,end)];
    s(2,:) = [res_before.H_ee(end-params.ne+1,end) res_before.H_ie(end-params.ni+1,end)];
    s(3,:) = [res_before.H_ei(end-params.ne+1,end) res_before.H_ii(end-params.ni+1,end)];
    s(4,:) = res_before.ref;
else
    s=zeros(4,params.ne+params.ni);
%     s(2,:) = [100*ones(1,params.ne) ones(1,params.ni)*50];
end
% s(1,:)=[ones(1,params.ne) 40*ones(1,params.ni)];
%s is state matrix. The first three rows indicate the triplet(V,H_e,H_i)for
%each neuron, the fourth row is 1 and 0 indicating whether a neuron is at refractory
%state(1 means at refractory).

m=zeros(3,params.ne+params.ni);
m(1,:)=c(1,:);
%m is mean matrix of the next random trial. The first row of m is the same as
%the first row of c, while the other three rows is obtained by c(2:4,:)./s(2:4,:)

res.spike=zeros(2000,params.ne+params.ni);
res.t1_Exsti = zeros(2000,params.ne+params.ni);
res.t2_Esti = zeros(2000,params.ne+params.ni);
res.t3_Isti = zeros(2000,params.ne+params.ni);
res.deltaV_Isti = zeros(2000,params.ne+params.ni);
%spike is the spike time train of each neuron. The first row of spike is spike count.


%% run model

while t<= duration_time
    switch need_set 
        case{'constant'} % constant possion
%             if t < 1000
%                 lambda = 2;
%             elseif t<2000
%                 lambda = 4;
%             elseif t<3000
%                 lambda = 6;
%             elseif t<4000
%                 lambda = 8;
%             elseif t<5000
%                 lambda = 10;
%             end
                
                
            lambda = params.Ex_Poisson_lambda;
            res.Ex_type = 'constant';
            
        case{'sin'} % sin possion
            
            C1 = params.Ex_sin_C1;
            C2 = params.Ex_sin_C2;
            C3 = params.Ex_sin_C3;
            C4 = params.Ex_sin_C4;
            if t < 1000
                C2 = 0.01*pi;
                lambda = C1*sin(C2*(t+C3))+C4;
            elseif t<2000
                C2 = 0.05*pi;
                lambda = C1*sin(C2*(t+C3))+C4;
            elseif t<3000
                C2 = 0.1*pi;
                lambda = C1*sin(C2*(t+C3))+C4;
            elseif t<4000
                C2 = 0.15*pi;
                lambda = C1*sin(C2*(t+C3))+C4;
            elseif t<5000
                C2 = 0.2*pi;
                lambda = C1*sin(C2*(t+C3))+C4;
            end

            res.Ex_type = 'sin';
    end

    Ex_input_tau = repmat(1/lambda, 1, params.ne + params.ni);
    m(1,:) = 1./(Ex_input_tau);
    m(2:4,:)=s(2:4,:)./c(2:4,:);
    sum_m=sum(sum(m));
    interval = exprnd(1/sum_m);
    p_array = cumsum(reshape((m/sum_m)',[],1));
    ind = sum(p_array<=rand(1))+1;
    x = floor((ind-1)/(params.ne+params.ni))+1;
    y = ind - (x-1)*(params.ne+params.ni);
    
    t_old = t;
    t = t + interval;
    
    % 增加采样时间点
    if fix(t_old/0.1) ~= fix(t/0.1)
        res.V_e = [res.V_e,s(1,1:params.ne)];
        res.V_i = [res.V_i,s(1,params.ne+1:params.ne+params.ni)];
        res.H_ee = [res.H_ee s(2,1:params.ne)];
        res.H_ei = [res.H_ei s(3,1:params.ne)];
        res.H_ie = [res.H_ie s(2,params.ne+1:params.ne+params.ni)];
        res.H_ii = [res.H_ii s(3,params.ne+1:params.ne+params.ni)];
        res.t = [res.t t];
        res.eref = [res.eref s(4,1:params.ne)];
        res.iref = [res.iref s(4,params.ne+1:params.ne+params.ni)];
    end
    
    
%     if fix(t_old*params.Ex_Poisson_lambda) ~= fix(t*params.Ex_Poisson_lambda)
%         for k = 1:params.ne+params.ni
%             deltaV = ceil((t*params.Ex_Poisson_lambda)-(t_old*params.Ex_Poisson_lambda));
%             res.t1_Exsti(1,k)=res.t1_Exsti(1,k)+deltaV;
%             res.t1_Exsti(res.t1_Exsti(1,k)+1,k)=t;            
%             if s(4,k)==0
%                 
% %                 res.t1_Exsti(1,k)=res.t1_Exsti(1,k)+deltaV;
% %                 res.t1_Exsti(res.t1_Exsti(1,k)+1,k)=t;
%                 s(1,k)=s(1,k)+deltaV;
%             end
%         end
%     end
   
    if x==1 % external input operates

        if s(4,y)==0
%             res.t1_Exsti(1,y)=res.t1_Exsti(1,y)+1; 
%             res.t1_Exsti(res.t1_Exsti(1,y)+1,y)=t;
            s(1,y)=s(1,y)+1;
        end


    elseif x==2 % H_e operates
        s(x,y)=s(x,y)-1;
%         res.t2_Esti(1,y)=res.t2_Esti(1,y)+1; 
%         res.t2_Esti(res.t2_Esti(1,y)+1,y)=t;
        if s(4,y)==0
            if y<= params.ne
                s(1,y)=s(1,y)+params.s_ee;
            else
                s(1,y)=s(1,y)+params.s_ie;
            end
        end
        
    elseif x==3 % H_i operates
%         res.t3_Isti(1,y)=res.t3_Isti(1,y)+1; 
%         res.t3_Isti(res.t3_Isti(1,y)+1,y)=t;
        res.deltaV_Isti(1,y)=res.deltaV_Isti(1,y)+1; 
        if y<=params.ne
            res.deltaV_Isti(res.deltaV_Isti(1,y)+1,y)=(s(1,y)+params.Mr)*params.s_ei/(params.M+params.Mr);
%             res.deltaV_Isti(res.deltaV_Isti(1,y)+1,y)=params.s_ei;
        else
            res.deltaV_Isti(res.deltaV_Isti(1,y)+1,y)=(s(1,y)+params.Mr)*params.s_ii/(params.M+params.Mr);
%             res.deltaV_Isti(res.deltaV_Isti(1,y)+1,y)=params.s_ii;
        end
        s(x,y)=s(x,y)-1;
        if s(4,y)==0
            if y<=params.ne
                sv=(s(1,y)+params.Mr)*params.s_ei/(params.M+params.Mr);
%                 sv=params.s_ei;
                ssv=floor(sv);
                if sv-ssv-rand(1)>0
                    s(1,y)=s(1,y)-ssv-1;
                else
                    s(1,y)=s(1,y)-ssv;
                end
            else
                sv=(s(1,y)+params.Mr)*params.s_ii/(params.M+params.Mr);
%                 sv=params.s_ii;
                ssv=floor(sv);
                if sv-ssv-rand(1)>0
                    s(1,y)=s(1,y)-ssv-1;
                else
                    s(1,y)=s(1,y)-ssv;
                end
            end
        end
    else % recover from refractory
        s(4,y)=0;
    end
    

    if s(1,y)>=params.M %a neuron reaches threshold and spike
        s(1,y)=0;
        s(4,y)=1;
        if y<=params.ne %an excitatory spike
            %decide the postsynaptic neurons
            p=rand(1,params.ne+params.ni);
            p(1:params.ne)=-p(1:params.ne)+params.p_ee; 
            p(params.ne+1:params.ne+params.ni)=-p(params.ne+1:params.ne+params.ni)+params.p_ie;
            p=sign(abs(p)+p);
            s(2,:)=s(2,:)+p;
        else %an inhibitory spike
            %decide the postsynaptic neurons
            p=rand(1,params.ne+params.ni);
            p(1:params.ne)=-p(1:params.ne)+params.p_ei; 
            p(params.ne+1:params.ne+params.ni)=-p(params.ne+1:params.ne+params.ni)+params.p_ii;
            p=sign(abs(p)+p);
            s(3,:)=s(3,:)+p;
        end
        res.spike(1,y)=res.spike(1,y)+1; %write down the spike time
        res.spike(res.spike(1,y)+1,y)=t;

   
    end
end

end
