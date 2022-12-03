function plot_fft(X1,X2,single)
if nargin > 3
    error('输入变量过多！');
elseif nargin == 2
    single = false; % 默认情况下为all
end
bin = 0.1;
duration_time = 1000;
%% spike数目统计图
% X1(1) = [];
% X2(1) = [];
% figure;
% plot(bin:bin:duration_time,X1)
% title('X')
% xlabel('t (milliseconds)')
% ylabel('X(t)')


%% rasterplot
ne = 300;
ni = 100;
spike = zeros(20000,ne+ni);
X1(X1<0) = 0;
X2(X2<0) = 0;
figure
for i = 1:duration_time/bin
    spike(i,randi(ne,1,X1(i)))=i*0.1;
    spike(i,ne+randi(ni,1,X2(i)))=i*0.1;
end 

for i=1:ne
times = spike(:,i);
times = times(times~=0);
num   = size(times, 1);
scatter(times, i*ones(num, 1),12,'.','r');
hold on
end

for i=(ne+1):(ne+ni)
times = spike(:,i);
times = times(times~=0);
num  = size(times, 1);
scatter(times, i*ones(num, 1),12,'.','b');
hold on
end

%% 频谱图
% M = reshape(X1,1/bin,[]);
% X1 = sum(M);
% L = duration_time;
% Y = fft(X1);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(1) = 0;
% P1(2:end-1) = 2*P1(2:end-1);
% f = (1000)*(0:(L/2))/L;
% [P1_sort, pos] = sort(P1, 'descend');
% fprintf('The frequencies of the top five components are %2.1f %2.1f %2.1f %2.1f %2.1f Hz, the amplitudes are %2.2f %2.2f %2.2f %2.2f %2.2f \n', f(pos(1:5)), P1_sort(1:5));
% 
% figure;
% plot(f,P1)
% title('Single-Sided Amplitude Spectrum of Spike Density')
% xlabel('f (Hz)')
% ylabel('amplitude')
end