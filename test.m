clear all;
close all;
clc;
%
% 参数初始化
fc=10000;   %载波频率
% fs=fc*3;  %采样频率
%fs=240000;    %接收端采样频率
fs=96000;%发送端采样频率，后面用接收器时要切换！！
fb=200;    %符号速率
% fb=560;    %符号速率
df=0;       %初始频偏
fi=fc+df;   %输入信号载波频率
% SNR=30;   %输入数据信噪比(dB)
ts=1/fs;
OverSamp=fs/fb;
Delay=5;
back=[1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];%13位巴克码
bitlen=100;
rand('seed',0);
x1 = randi([0,1],1,bitlen);
x1 = sign(x1-0.5); % 01码元转换为正负一
x2 = randi([0,1],1,bitlen);
x2 = sign(x2-0.5);
x3 = randi([0,1],1,bitlen);
x3 = sign(x3-0.5);
s_send = [back back back x3]; % 发送序列
plot(s_send);

t = readtable('bpsk_10k_96k.dat', 'HeaderLines', 1);
a_0 = table2array(t);
threshold = max(abs(a_0)) * 0.5; % 设定阈值为最大幅度的一半，这样非常稳健，不怕底噪
% 找到最后一个大于阈值的点（即高电平的结束点）
last_high_idx = find(abs(a_0) > threshold, 1, 'last');
% 截取数据：只保留到高电平结束的那一点
a = a_0(1:last_high_idx);
 %a = a(1:min(100000, length(a)), :); 
a = a(1:length(a)-1,:);%无条件地删除信号 a 的最后一个数据点（样本）。
figure;plot(a)
title('原始接收信号波形');
xlabel('样本点'); ylabel('幅度');
Signal = a;
nn=length(Signal);
f_s=(-nn/2:nn/2-1)*(fs/nn);
figure;plot(f_s,abs(fftshift(fft(Signal))))
title('原始接收信号频谱');
xlabel('频率 (Hz)'); ylabel('幅度');

if isrow(Signal)%强转成列向量
    Signal = Signal(:); 
end
%% 二阶带通滤波器
% 设定带通滤波器的频率范围
f_low = 9000;   % 下限频率为9.98KHz
f_high = 11000; % 上限频率为10.2KHz（确保信号包含1.02KHz成分）

% 计算滤波器的归一化频率
Wn = [f_low f_high] / (fs / 2);  % MATLAB中频率归一化范围是[0,1]，fs/2是奈奎斯特频率

% 使用butter设计带通滤波器（4阶滤波器）
[b, a] = butter(4, Wn, 'bandpass');

% 对信号进行带通滤波双向滤波，相位延迟为0
y_filtered = filtfilt(b, a, Signal);

% 进行滤波后的信号频域显示
NLen = length(y_filtered);  % 信号长度
f_s = (-NLen/2 : NLen/2-1) * (fs / NLen);  % 频率范围

% 频谱计算
Y_filtered_fft = fft(y_filtered, NLen);
Y_filtered_mag = abs(fftshift(Y_filtered_fft));  % 对频谱进行频移以便显示

% 绘制带通滤波后信号的频谱
figure;
plot(f_s, Y_filtered_mag);
title('带通滤波后的频谱');
xlabel('频率 (Hz)');
ylabel('幅度');

% 绘制带通滤波后的时域信号
figure;
plot(y_filtered);
title('带通滤波后的时域信号');
xlabel('采样点');
ylabel('幅度');

num=length(y_filtered)-1;  %数据长度
o_sin=sin(2*pi*fc*(0:num)*ts); % 本地同相参考
din=y_filtered;
alpha=0.8;
Delay=5;
din = din(:); 
o_sin=o_sin(:);
% 解调（假设 BPSK：用 I 路即可）
DemodWave=din.*o_sin;                                 % 用同步后的本地正弦乘，得到基带 BPSK（同相路）
plot(f_s,abs(fftshift(fft(DemodWave))))  ;
h_sqrt=rcosine(1,OverSamp,'fir/sqrt',alpha,Delay);    % 根升余弦（RRC）滤波器（旧函数；新版本可用 rcosdesign）
RcvMatched=conv(DemodWave,conj(h_sqrt));              % 匹配滤波（卷积实现；也可用 filter 并考虑群时延补偿）

%%
% 频谱与时域可视化
figure;
nn=length(RcvMatched);
f_s=(-nn/2:nn/2-1)*(fs/nn);                            % 频率轴（Hz）
plot(f_s,abs(fftshift(fft(RcvMatched))))               % 匹配滤波后信号幅度谱
figure;plot(DemodWave);title('解调后的时域波形')
figure;plot(RcvMatched);title('匹配滤波后的时域波形')


%%
%第一次降采样（过采样率是480）
D1=120;
% 这里的 '1' 是假设的起始点。在实际工程中，可能需要根据滤波器延迟调整这个起始点
Signal_4sps = RcvMatched(1 : D1 : end);
figure;
plot(Signal_4sps,'.-');
title(['降低过采样率后,原采样率降低倍数: ' num2str(D1)]); 

%%
%载波同步




%%
%二次降采样,此时要找最佳采样点
D2=4;
% 我们不能盲目地取第1个，或者求平均。我们要找能量最大的那个相位。
max_pwr = 0;
best_phase = 1;
for phase = 1 : D2
    % 每隔4个点取一个，尝试不同的起始位置
    temp_signal = Signal_4sps(phase : D2 : end);
    pwr = mean(abs(temp_signal).^2); % 计算能量
    if pwr > max_pwr
        max_pwr = pwr;
        best_phase = phase;
    end
end
fprintf('自动检测到的最佳采样相位是: %d\n', best_phase);
% 执行最终抽取
Signal_1sps = Signal_4sps(best_phase : D2 : end);
figure;
stem(Signal_1sps);
title(['最终符号 (1 sps) - 最佳相位: ' num2str(best_phase)]);
%%
%位同步


%%
%计算相位星座图
plot(real(Signal_1sps), zeros(size(Signal_1sps)), '.');
title('BPSK 星座图 (实数)');
xlabel('In-Phase (I)'); ylabel('Quadrature (Q)');
axis([-2 2 -1 1]); % 限制范围方便观察

%%
%符号判决
Soft_Bits = real(Signal_1sps);
% 映射规则要与发送端一致。假设发送端 1->1, 0->-1
RcvBit = zeros(1, length(Soft_Bits));
for k = 1:length(Soft_Bits)
    if Soft_Bits(k) > 0
        RcvBit(k) = 1;   % 正电平判为 1
    else
        RcvBit(k) = -1;  % 负电平判为 -1
    end
end
%% 误比特率
%% conv算法
selfc_send = conv(s_send,fliplr([back back back]));% fliplr用于翻转向量
[a_send ,index_send]=max(abs(selfc_send))
figure;stem(selfc_send); title('发送序列与本地PN序列的互相关函数');
pian_send = abs(length([back back back])-index_send)
s_recive=RcvBit;
selfc_rece = conv(s_recive,fliplr([back back back]));% fliplr用于翻转向量
[a_rece ,index_rece]=max(abs(selfc_rece))
figure;stem(selfc_rece); title('接收序列与本地PN序列的互相关函数');
pian_rece = abs(length([back back back])-index_rece)
%%
% 误比特率 
s_send_ch=sign(s_send+1);
s_recive_ch=sign(s_recive+1);
len=bitlen;% 计算800个bit的误比特率
% 计算可访问的最大长度
max_len = min(length(s_send_ch) - pian_send, length(s_recive_ch) - pian_rece);

if len > max_len
    warning('bitlen超出有效比特长度，自动调整为最大可用长度 %d', max_len);
     len = max_len-1;
end
x=s_send_ch(pian_send+1:pian_send+1+len);
y=s_recive_ch(pian_rece+1:pian_rece+1+len);
x = x(:); % 强制转换为列向量
y = y(:); % 强制转换为列向量
x(x == -1) = 0;
[number,ratio] = biterr(x,y)% number误比特的个数，ratio：误比特率
%% matlab互相关算法
% 发送序列与巴克码的互相关函数

s_recive=a;
[selfc_send,track_send]=xcorr(s_send,[back back back]);% self_c:互相关函数,track:延迟矢量，用于画图
figure(11); 
selfc_send=abs(selfc_send);
plot(selfc_send)
title('发送序列自相关函数');
%
[selfc_rece,track_rece]=xcorr(s_recive,[back back back]);% self_c:互相关函数,track:延迟矢量，用于画图
figure(12); 
selfc_rece=abs(selfc_rece);
plot(selfc_rece)
title('接收序列与本地序列的互相关函数');
%%
% 求同步偏移量
[a_send ,index_send]=max(selfc_send);
[a_rece ,index_rece]=max(selfc_rece);
pian_rece = abs(length(s_recive)-index_rece);
pian_send = abs(length(s_send)-index_send);



%% 误码分布图
E_N=zeros(1,length(x));
Err_num=0;
for i=1:length(x)
    if x(i)~=y(i)
        E_N(i)=1;
        Err_num=Err_num+1;
    end
end
figure;plot(E_N)


