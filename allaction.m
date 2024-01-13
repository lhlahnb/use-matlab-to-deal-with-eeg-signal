clear all 
clc
data = load('EEG_Data.mat');
s5_bis40 = data.s5_bis40;
s5_bis60 = data.s5_bis60;
s5_wake = data.s5_wake;
s = [s5_bis40;s5_bis60;s5_wake];
length = 2*256*60;
% 任务一 使用pwelch进行EEG信号的频谱分析
%按照ppt设置参数
NFFT= 1024;%分段数据点数
Fs = 256;%数据采样率
Window = hanning(1250);%窗的类型：汉宁窗
Noverlap = 0.9;%数据重叠率
%进行处理
[Pxx1,fx1] = pwelch(s(1,:)',Window,Noverlap,NFFT,Fs,'psd');
[Pxx2,fx2] = pwelch(s(2,:)',Window,Noverlap,NFFT,Fs,'psd');
[Pxx3,fx3] = pwelch(s(3,:)',Window,Noverlap,NFFT,Fs,'psd');
%绘图
figure
subplot(3,1,1)
plot(fx1,Pxx1)
xlabel('Frequency/Hz')
title("power spectrum density--bis40");
subplot(3,1,2)
plot(fx2,Pxx2)
xlabel('Frequency/Hz')
title("power spectrum density--bis60");
subplot(3,1,3)
plot(fx3,Pxx3)
xlabel('Frequency/Hz')
title("power spectrum density--wake");
% 任务二 EEG信号各频段的信号提取；
%设计一个五阶的butterworth滤波器
%通带边界为 1-4,4-8,8-13,13-30,30-45 单位Hz
filtkernel = cell(1,5);
filtkernel{1} = designfilt('bandpassiir','FilterOrder',18, ...
'HalfPowerFrequency1',1,'HalfPowerFrequency2',4, ...
'SampleRate',256);
filtkernel{2} = designfilt('bandpassiir','FilterOrder',18, ...
'HalfPowerFrequency1',4,'HalfPowerFrequency2',8, ...
'SampleRate',256);
filtkernel{3} = designfilt('bandpassiir','FilterOrder',18, ...
'HalfPowerFrequency1',8,'HalfPowerFrequency2',13, ...
'SampleRate',256);
filtkernel{4} = designfilt('bandpassiir','FilterOrder',18, ...
'HalfPowerFrequency1',13,'HalfPowerFrequency2',30, ...
'SampleRate',256);
filtkernel{5} = designfilt('bandpassiir','FilterOrder',18, ...
'HalfPowerFrequency1',30,'HalfPowerFrequency2',45, ...
'SampleRate',256);
%只使用清醒和BIS40的数据
%进行滤波
s5_bis40_band = zeros(5,length);
s5_wake_band = zeros(5,length);
for i = 1:5
s5_bis40_band(i,:) = filter(filtkernel{i},s(1,:));
s5_wake_band(i,:) = filter(filtkernel{i},s(3,:));
end
%绘图 前五秒
l = 5*Fs;
t = linspace(0,5,l);
figure
strs = {'s4-bis40','s4-bis60','s4-wake'};
strband = {'Delta','Theta','Alpha','Beta','Gamma'};
subplot(5,2,1)
plot(t,s5_bis40_band(1,1:l))
title("Delta-BIS40");
xlabel('Time/s')
subplot(5,2,2)
plot(t,s5_wake_band(1,1:l))
title("Delta-wake")
xlabel('Time/s')
subplot(5,2,3)
plot(t,s5_bis40_band(2,1:l))
title("Theta-BIS40");
xlabel('Time/s')
subplot(5,2,4)
plot(t,s5_wake_band(2,1:l))
title("Theta-wake")
xlabel('Time/s')
subplot(5,2,5)
plot(t,s5_bis40_band(3,1:l))
title("Alpha-BIS40");
xlabel('Time/s')
subplot(5,2,6)
plot(t,s5_wake_band(3,1:l))
title("Alpha-wake")
xlabel('Time/s')
subplot(5,2,7)
plot(t,s5_bis40_band(1,1:l))
title("Beta-BIS40");
xlabel('Time/s')
subplot(5,2,8)
plot(t,s5_wake_band(1,1:l))
title("Beta-wake")
xlabel('Time/s')
subplot(5,2,9)
plot(t,s5_bis40_band(1,1:l))
title("Delta-BIS40");
xlabel('Time/s')
subplot(5,2,10)
plot(t,s5_wake_band(1,1:l))
title("Delta-wake")
xlabel('Time/s')
% 任务三 EEG信号各波段能量占比
%计算功率比
power_pwelch = zeros(3,6);
power_pwelch(1,:)=powerRate(Pxx1,fx1); 
power_pwelch(2,:)=powerRate(Pxx2,fx2); 
power_pwelch(3,:)=powerRate(Pxx3,fx3); 
%绘图
figure
Xstr = categorical({'Slow','Delta','Theta','Alpha','Beta','Gamma'});
Xstr = reordercats(Xstr,{'Slow','Delta','Theta','Alpha','Beta','Gamma'});
subplot(3,1,1)
bar(Xstr,power_pwelch(1,:))
title("power ratio--BIS40")
xlabel('Frequency/Hz');
ylabel('Percent %');
subplot(3,1,2)
bar(Xstr,power_pwelch(2,:))
title("power ratio--BIS60")
xlabel('Frequency/Hz');
ylabel('Percent %');
subplot(3,1,3)
bar(Xstr,power_pwelch(3,:))
title("power ratio--wake")
xlabel('Frequency/Hz');
ylabel('Percent %');
% 任务四 随麻醉状态变化的EEG信号各波段占比
data=[s5_wake,s5_bis60,s5_bis40]; 
NFFT=1024; 
Fs=256; 
Window=hanning(1250); 
Noverlap=0.9; 
N_win=size(data,2)/(5*Fs); 
power=zeros(6,N_win); 
for i=1:N_win 
[Pxx,fx]=pwelch(data(1+5*Fs*(i-1):5*Fs+5*Fs*(i-1))',Window,Noverlap,NFFT,Fs,'power'); 
power(:,i)=powerRate(Pxx,fx); 
end
% 绘图
figure
for i = 1:6
plot(1:N_win,power(i,:),'-o')
hold on 
end
xlabel("Window Number")
ylabel("Percent%");
title("power ratio-");
legend('Slow','Delta','Theta','Alpha','Beta','Gamma');
% 任务五 EEG信号的时频分析
tas5 = spectrogram(data);
figure
wind = hanning(256);
nfft = 256;

spectrogram(data,wind,1,nfft,Fs,'yaxis')


%% powerRate函数
function [power]=powerRate(Pxx,fx) 
all_idx=find(fx>=0.1 & fx<=45);
slow_idx=find(fx>=0.1 & fx<=1); 
slow=100*sum(Pxx(slow_idx))/sum(Pxx(all_idx)); 
delta_idx=find(fx>1 & fx<=4); 
delta=100*sum(Pxx(delta_idx))/sum(Pxx(all_idx)); 
theta_idx=find(fx>4 & fx<=8); 
theta=100*sum(Pxx(theta_idx))/sum(Pxx(all_idx)); 
alpha_idx=find(fx>8 & fx<=13); 
alpha=100*sum(Pxx(alpha_idx))/sum(Pxx(all_idx)); 
beta_idx=find(fx>13 & fx<=30); 
beta=100*sum(Pxx(beta_idx))/sum(Pxx(all_idx)); 
gamma_idx=find(fx>30 & fx<=45); 
gamma=100*sum(Pxx(gamma_idx))/sum(Pxx(all_idx)); 
power=[slow;delta;theta;alpha;beta;gamma]; 
end