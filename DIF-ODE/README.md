# 命令行
## LIF模型：
\>main_LIF
\>plot(params.dt:params.dt:params.duration_time,res_lif.fr_e,'r');
## DIFODE模型：
## 无初始状态，havestart设置为0
\>main_DIFODE     
\>plot(params.dt:params.dt:params.duration_time,res_DIFODE.fr_e,'r');
## 有初始状态，havestart设置为1
\>get_start     
\>main_DIFODE     
## 绘制各个曲线神经元数目的动态变化动画
\>plot_ODEn
