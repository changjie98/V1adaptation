# 命令行
## LIF模型：
\>main_LIF     
\>plot(params.dt:params.dt:params.duration_time,res_lif.fr_e,'r');     % 观察放电率变化
\>plot(params.dt:params.dt:params.duration_time,mean(res_lif.V_e,2),'r');     % 观察平均电压变化

## DIFODE模型：
## 无初始状态，havestart设置为0
\>main_DIFODE     
\>plot(params.dt:params.dt:params.duration_time,res_DIFODE.fr_e,'r');     
\>plot(params.dt:params.dt:params.duration_time,mean(res_DIFODE.V_e,2),'r');     % 观察平均电压变化
## 有初始状态，havestart设置为1
\>get_start     
\>main_DIFODE     
## 绘制各个曲线神经元数目的动态变化动画
\>plot_ODEn
