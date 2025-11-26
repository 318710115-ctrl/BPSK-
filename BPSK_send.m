clear;
close all;
clc;
%%
% 参数初始化
fc=14000;%载波频率
fb=875;%符号速率
fs=336000;%采样频率 
% fs=fc*4*6;%采样频率 
OverSamp=fs/fb;%过采样率
%% 生成发送码元
back=[1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];%13位巴克码
bitlen=800;
rand('seed',0);
x=randi([0,1],1,bitlen);
x=sign(x-0.5);%01码元转换为正负一
% x=[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,...
%   -1,-1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,1,-1,-1,1,-1,1,-1,-1,1,-1,1,-1,1,1,-1,1,1,1,...
%   1,-1,1,1,1,1,1,-1,-1,-1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,1,1,1,-1,1,1,-1,1,-1,-1,1,...
%   1,1,-1,-1,-1,-1,1,-1,-1,1];
s_send=[back back back x];
%%
% 脉冲成型
% A,B 卷积后的长度为length(A)+length(B)-1
Delay=5;
alpha=0.8;  
h_sqrt=rcosine(1,OverSamp,'fir/sqrt',alpha,Delay);% h_sqrt为滤波器的冲激响应，长度=Delay*OverSamp*2+1
SendSignal_OverSample=kron(s_send,[1 zeros(1,OverSamp-1)]);% 过采样的发送信号
SendShaped=conv(SendSignal_OverSample,h_sqrt);   % 余弦滚降滤波后的采样信号
%% 基带信号分析
%波形
figure;plot(SendShaped);
%频谱
nn=length(SendShaped);
f_s=(-nn/2:nn/2-1)*(fs/nn);
figure;plot(f_s,abs(fftshift(fft(SendShaped))))
%% 调制
N=0:length(SendShaped)-1;
CarrierWave=sin(2*pi*fc*N/fs); % 载波波形
ModemWave=round((SendShaped.*CarrierWave).*100000000);% 频带信号
figure;plot(ModemWave)
% 频带信号频谱
nn=length(ModemWave);
f_s=(-nn/2:nn/2-1)*(fs/nn);
figure;plot(f_s,abs(fftshift(fft(ModemWave))))
% x=textread('E:\Mission\mission12\a5.txt','%s');
% for i = 1:length(x)
%     a(i) = str2num(x{i});
% end
% ModemWave=a;
%%
save 2psk_fc_14_fb875_326016 ModemWave
dlmwrite('2p.txt',ModemWave,'delimiter',' ','precision','%2.4f')%空格间隔
%%
% 发声
% sound(ModemWave,fs);
%% 接收端
Signal=ModemWave;
NSym=ceil(length(Signal)/OverSamp);
alpha=0.8;
df=0;    %初始频偏
fi=fc+df;   %输入信号载波频率
ts=1/fs;
K=4;     %每个符号采4个样点
Ns=K*NSym;  %总的采样点数
%%
% 50Hz陷波器 
f0=50;   %工频干扰频率为50HZ 
Ts=ts;   %采样间隔为0.001 
%fs=1/Ts;   
NLen=length(Signal);    %长度为1602000
n=0:NLen-1;  %陷波器的设计 
apha=-2*cos(2*pi*f0*Ts); 
beta=0.96; 
b=[1 apha 1]; 
a=[1 apha*beta beta^2]; 
% figure(200); 
% freqz(b,a,NLen,fs);          %陷波器特性显示 
x=Signal;  %原始信号 
y0=dlsim(b,a,x);                     %陷波器滤波处理 
xfft=fft(x,NLen); xfft=xfft.*conj(xfft)/NLen;         %对信号进行频域变换 
y1=fft(y0,NLen); y2=y1.*conj(y1)/NLen; 

figure
% plot(y0(2000:9000)');
nn=length(y0');
f_s=(-nn/2:nn/2-1)*(fs/nn);
figure(1);plot(f_s,abs(fftshift(fft(y0'))))
title('50Hz陷波后的信号频谱');

%% 二阶带通滤波器
% 设定带通滤波器的频率范围
f_low = 9000;   % 下限频率为9.98KHz
f_high = 11000; % 上限频率为10.2KHz（确保信号包含1.02KHz成分）

% 计算滤波器的归一化频率
Wn = [f_low f_high] / (fs / 2);  % MATLAB中频率归一化范围是[0,1]，fs/2是奈奎斯特频率

% 使用butter设计带通滤波器（2阶滤波器）
[b, a] = butter(2, Wn, 'bandpass');

% 对信号进行带通滤波
y_filtered = filter(b, a, Signal);

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

%%
% 载波同步
num=length(Signal)-1;  %数据长度
%本地VCO的增益，也相当于最大捕获频偏
K0=10e3; %VCO增益，即最大捕获频偏   
% 载波同步环路滤波器系数计算
BL=0.05*fb;
el=0.707;
% K=1373.75;
T=10/fs;
N_N=6;                        %数据位宽
Bip=6;
% K=(2*pi*fs)/(2^N_N)*T*(2^(Bip-2));
K_zbtb=1;
wn=BL*((8*el)/(1+4*el));
wn1=wn/(2*pi);
c1=((2*el*wn1*T)/K_zbtb);          %环路滤波器系数c1
c2=((wn1*T)^2/K_zbtb);             %环路滤波器系数c2
% %环路滤波器系数
% c1=2^(-14);
% c2=2^(-19);
%环路输入信号
%din=Signal ;
din=y0';
%%
% 设计 IIR 低通滤波器，用 Costas 环处理信号
N=5;                 % IIR 低通滤波器阶数（Chebyshev II 型）
R=60;                % 阻带衰减 (dB)
Wn=15000;            % （模拟）截止频率，单位 Hz
[lpf_b,lpf_a]=cheby2(N,R,Wn*2/fs);   % 设计 Chebyshev-II IIR 低通；Wn 需归一化到 [0,1]（相对 fs/2），因此用 Wn*2/fs

% 估计滤波器的增益（示例：取分子/分母首项之比作为 DC 增益近似）
% glpf=sum(abs(lpf_b))            % 另一种粗略写法：零点系数绝对值之和
glpf=lpf_b(1)/lpf_a(1);          % 这里取 b(1)/a(1) 作为近似

% 本地载波（NCO）输出信号初值（先按固定 fc 生成的参考正弦/余弦）
o_sin=sin(2*pi*fc*(0:num-1)*ts); % 本地同相参考
o_cos=cos(2*pi*fc*(0:num-1)*ts); % 本地正交参考

% 乘法器（下变频）输出缓存初始化
mult_i=zeros(1,num);              % I 路（与 sin 相乘）
mult_q=zeros(1,num);              % Q 路（与 cos 相乘）

% 鉴相器（相位误差）输出缓存初始化
pd=zeros(1,num);                  % 存储每个采样点的相位误差估计

% 本地载波瞬时频率轨迹初始化（起始为 fc）
fo=ones(1,num)*fc;

% —— 开始 Costas 载波同步环路处理 —— %
% 对 Costas 环路使用到的状态变量初始化
len=0;                 % 低通输出缓冲长度（用于取最后一个样本）
dfreq=0;               % 环路滤波器输出（频率/相位控制量，馈入 VCO）
temp=0;                % 环路滤波器的积分状态（PI 的积分项）
df=0;                  % 每次更新时 VCO 的输入（相当于频偏/相位误差的控制量）
thera=0;               % 本地 VCO 输出信号的“相位偏移”（用于块间相位连续）
n=2;                   % 记录 VCO 更新次数（用于 thera 索引）
lvco=10;               % VCO 频率的更新周期：每 lvco 个采样点更新一次 df
m=0;                   % 周期内的样本计数（0~lvco-1）

for i=(lvco+1):(num)   % 从第 lvco+1 个样本开始处理（保证有足够的历史样本用于滤波）
   
    if (din(i)==0)
        din(i)=10e-30; % 避免完全为 0 导致后续除零或 atan 不稳定（加入极小量）
    end
    
    % 下变频：分别与本地 sin/cos 相乘，得到 I/Q 两路
    mult_i(i)=din(i).*o_sin(i);
    mult_q(i)=din(i).*o_cos(i);
    
    % 低通滤波：提取基带分量
    % 注意：这里每次仅对一小段 [i-N-2 : i] 进行 filter，相当于频繁重置滤波器状态，
    % 在严格意义上不等价于“持续的 IIR”，但计算量小；若要精确，应维护滤波器状态 zi/zo。
    lpf_i=filter(lpf_b,lpf_a,mult_i(i-N-2:i));
    lpf_q=filter(lpf_b,lpf_a,mult_q(i-N-2:i));
   
    % 鉴相：以 I/Q 低通后“最新值”计算相位误差
    len=length(lpf_i);                 % 当前短窗长度
    pd(i)=atan(lpf_q(len)/lpf_i(len)); % 误差 = atan(Q/I)；更稳健可用 atan2(Q,I)

    % 环路滤波器（PI 结构）：dfreq = Kp*e + 累积积分
    dfreq = c1*pd(i)+temp;             % 比例 + 积分状态 = “控制量”候选
    temp = temp+c2*pd(i);              % 更新积分状态（累加误差）

    % 每 lvco 个采样点更新一次 VCO 的驱动量（降频更新以省计算）
    if(mod(i,lvco)==0)
        df=dfreq;                      % 锁存当前滤波器输出，作为接下来一个块的 df（频率偏差）
        % 计算从上次更新到本次更新这段的“相位累计”，以保证相位连续
        % thera(n-1) 是上一个块的相位偏移；2*pi*lvco*ts*fo(i-1) 是本块频率在 lvco 点内的相位累计
        thera(n)=2*pi*lvco*ts*fo(i-1)+thera(n-1);
        n=n+1;                         % 块计数 +1
    end

    % 更新本地载波瞬时频率（以 df 为偏移；K0 为 VCO 增益）
    fo(i)=fc+K0*df;

    % 计算当前样本在块内的位置（0~lvco-1）
    m=mod(i,lvco);

    % 生成“相位连续”的本地 VCO 输出载波（下一个时刻用）
    % 相位 = 当前频率 * 块内样本时间 + 上一次块边界相位偏移
    % 注意：这里写 o_sin(i+1)/o_cos(i+1) 可能在末尾 i=num 时越界；实际使用需注意边界保护
    o_sin(i+1)=sin(2*pi*fo(i)*m*ts+thera(n-1));
    o_cos(i+1)=cos(2*pi*fo(i)*m*ts+thera(n-1));
end

% 绘制：输入信号真实载波频率 vs NCO 跟踪到的频率
t=(0:num-1)*ts*1000;       % 时间轴（毫秒）
tfi=fi*ones(1,num);        % 输入信号的真实频率轨迹（常数）
figure(4);
plot(t,tfi/1e6,'r-.',t,fo/1e6,'b-');  % 红：真实；蓝：NCO 输出
legend('Frequency of input signal','Frequency of NCO output signal');
xlabel('Time(ms)');
ylabel('Frequency(MHz)');
grid on;

%%
% 解调（假设 BPSK：用 I 路即可）
DemodWave=din.*o_sin;                                 % 用同步后的本地正弦乘，得到基带 BPSK（同相路）
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
% —— 定时恢复与采样判决（Gardner + Farrow）——

% 环路滤波器输出（控制 NCO 步长/相位的量）寄存器；初值 0.5
w=[0.5,zeros(1,NSym-1)];

% NCO 相位寄存器（分数间隔累加器），初值 0.7（表示当前分数位置）
nco=[0.7 zeros(1,Ns-1)];

% NCO 临时状态（保存每步的 nco 减去步长后的中间值）
n_temp=[nco(1),zeros(1,Ns-1)];

% 分数间隔序列 u（0~1）：插值器所需的分数相位；初值 0.6
u=[0.6,zeros(1,2*NSym-1)];

% 插值输出缓存：I/Q 两路
yI=zeros(1,2*NSym);
yQ=zeros(1,2*NSym);

% Gardner 提取的定时误差缓存
time_error=zeros(1,NSym);

% 各类时间索引（离散“计数器”）
Ts_i=1;    % 输入采样的时间序号（用于索引 aI/bQ、n/n_temp/nco）
Ts_k=1;    % 插值输出的时间序号（用于索引 u/yI/yQ）
ms=1;      % 符号时序计数（用于 time_error 和 w 的更新）

strobe=zeros(1,Ns);                                  % “符号抽头”指示：用于判断何时计算 Gardner 误差

% 环路滤波器 PI 系数（可调，影响收敛和噪声）
% C11=5.41e-1; C22=3.82e-3;                           % 一组较大的带宽配置（注释掉）
C11=1.7e-3;   C22=1.5e-6;                              % 一组较小带宽配置（更平滑，收敛较慢）

%
%%%%% 仿真输入测试的 PSK 基带数据（形成插值输入序列） %%%%%
Len=length(RcvMatched);                                % 接收数据长度

% 生成插值输入序列（一般从匹配滤波输出中抽取过采样点）
cc=OverSamp;
rate=OverSamp/cc;                                      % 这里通常 = 1
aI=[RcvMatched(0*rate+1:rate:Len),0,0];                % I 路输入序列（末尾补 0 防越界）
bQ=[RcvMatched(0*rate+1:rate:Len),0,0];                % Q 路输入序列

% 可选：bQ=bQ1;                                        % 如果有单独的 Q 序列，这里可替换

% 主循环上界（-2 为了后续访问 Ts_i+2 不越界）
ns=length(aI)-2;

while(Ts_i<ns)

    % —— NCO 相位推进：用上一次的 nco 与当前步长 w(ms) 做差，得到下一时刻的临时相位 —— %
    n_temp(Ts_i+1)=nco(Ts_i)-w(ms);

    if(n_temp(Ts_i+1)>0)
        % 仍未跨过“符号边界”（相位未归零）：仅把临时值保存回 nco，继续累加
        nco(Ts_i+1)=n_temp(Ts_i+1);

    else
        % 已跨过“符号边界”（应当输出一个“符号采样点”）：
        % 1) 将 nco 相位取模，回到 [0,1)；2) 进行分数间隔插值，输出 yI/yQ；3)（必要时）计算 Gardner 误差并更新环路
        nco(Ts_i+1)=mod(n_temp(Ts_i+1),1);

        % —— 内插滤波器模块（Farrow 三次多项式形式）——
        % 基于 aI 的 4 点邻域（Ts_i-1, Ts_i, Ts_i+1, Ts_i+2）构造多项式系数：
        FI1=0.5*aI(Ts_i+2)-0.5*aI(Ts_i+1)-0.5*aI(Ts_i)+0.5*aI(Ts_i-1);
        FI2=1.5*aI(Ts_i+1)-0.5*aI(Ts_i+2)-0.5*aI(Ts_i)-0.5*aI(Ts_i-1);
        FI3=aI(Ts_i);
        % 根据当前分数间隔 u(Ts_k) 计算 I 路插值输出
        yI(Ts_k)=(FI1*u(Ts_k)+FI2)*u(Ts_k)+FI3;

        % Q 路同理
        FQ1=0.5*bQ(Ts_i+2)-0.5*bQ(Ts_i+1)-0.5*bQ(Ts_i)+0.5*bQ(Ts_i-1);
        FQ2=1.5*bQ(Ts_i+1)-0.5*bQ(Ts_i+2)-0.5*bQ(Ts_i)-0.5*bQ(Ts_i-1);
        FQ3=bQ(Ts_i);
        yQ(Ts_k)=(FQ1*u(Ts_k)+FQ2)*u(Ts_k)+FQ3;

        % 产生“抽头”标志：这里用 Ts_k 的奇偶来指示（每两个插值点当作一个符号）
        strobe(Ts_k)=mod(Ts_k,2);

        % —— 时钟误差提取（Gardner TED）——
        if(strobe(Ts_k)==0)
            % 每个“符号”计算一次误差（需要中心、中点、下一中心的样本）
            if(Ts_k>2)
               % 标准 Gardner：e = y_mid * (y_k - y_{k-1})（此处索引对应关系取决于具体实现）
               time_error(ms)=yI(Ts_k-1)*(yI(Ts_k)-yI(Ts_k-2)) + ...
                              yQ(Ts_k-1)*(yQ(Ts_k)-yQ(Ts_k-2));
            else
               % 前两个点不足以构造完整差分时，退化为一阶乘积
               time_error(ms)=(yI(Ts_k-1)*yI(Ts_k)+yQ(Ts_k-1)*yQ(Ts_k));
            end

            % —— 环路滤波器（PI）：用当前误差与差分更新步长控制量 w —— %
            if(ms>1)
                % w(k+1) = w(k) + Kp*(e(k)-e(k-1)) + Ki*e(k)
                w(ms+1)=w(ms)+C11*(time_error(ms)-time_error(ms-1)) + C22*time_error(ms);
            else
                % 第一拍没有 e(k-1)，等效只用 e(k)
                w(ms+1)=w(ms)+C11*time_error(ms) + C22*time_error(ms);
            end
            ms=ms+1;  % 符号计数 +1
        end

        % 插值输出计数 +1；并更新下一个插值点的分数间隔 u
        Ts_k=Ts_k+1;
        u(Ts_k)=nco(Ts_i)/w(ms);   % 分数相位比值（供下一次插值使用）
    end

    % 输入采样计数 +1
    Ts_i=Ts_i+1;
end
% y=yI+yQ;
  %%
  figure;
subplot(311);plot(u);xlabel('运算点数');ylabel('分数间隔');%axis([0 3.5e4 -1 1]);
subplot(312);plot(time_error);xlabel('运算点数');ylabel('定时误差');%axis([0 18e3 -0.5 0.5]);
subplot(313);plot(w);xlabel('运算点数');ylabel('环路滤波器输出');%axis([0 18e3 0.498 0.502]);
 %%
 % 符号判决
vv = OverSamp / 2;  % 原 cc/2
RcvBit = zeros(1, NSym);
max_i = floor(length(yI) / vv) * vv;  % 防止越界

for i = 1:vv:max_i
    if yI(i) > 0
        RcvBit((i - 1) / vv + 1) = -1;
    else
        RcvBit((i - 1) / vv + 1) = 1;
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
x(x == -1) = 0;
[number,ratio] = biterr(x,y)% number误比特的个数，ratio：误比特率
%% matlab互相关算法
% 发送序列与巴克码的互相关函数
s_recive=RcvBit;
[selfc_send,track_send]=xcorr(s_send,[back back back]);% self_c:互相关函数,track:延迟矢量，用于画图
figure(11); 
selfc_send=abs(selfc_send);
plot(selfc_send)
title('发送序列与接收序列的互相关函数');
%
[selfc_rece,track_rece]=xcorr(s_recive,[back back back]);% self_c:互相关函数,track:延迟矢量，用于画图
figure(12); 
selfc_rece=abs(selfc_rece);
plot(selfc_rece)
title('接收序列与接收序列的互相关函数');
%%
% 求同步偏移量
[a_send ,index_send]=max(selfc_send);
[a_rece ,index_rece]=max(selfc_rece);
pian_rece = abs(length(s_recive)-index_rece);
pian_send = abs(length(s_send)-index_send);
%% yI
yi_panjue=yI(1:vv:vv*(NSym));
figure;plot(yi_panjue,'b.')
figure;plot(yi_panjue(pian_rece:end),'b.')
%% RcvMatch
figure;plot(RcvMatched)
figure;plot(RcvMatched((pian_rece)*OverSamp:(pian_rece+39)*OverSamp))
%% DemodWave
figure;plot(DemodWave)
figure;plot(DemodWave((pian_rece-Delay)*OverSamp:(pian_rece+39-Delay)*OverSamp))
%% Signal
figure;plot(Signal)
figure;plot(Signal((pian_rece-Delay)*OverSamp:(pian_rece+39-Delay)*OverSamp))
%%
% %误码率曲线图
% for len=1:90000
%     x=s_send(pian_send+1:pian_send+1+len);
%     y=s_recive(pian_rece+1:pian_rece+1+len);
%     [number,ratio] = biterr(x,y);% number误比特的个数，ratio：误比特率
%     BER(len)=ratio;
% end
% figure(13);
% plot(BER);

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
%% 一些测试
% figure;plot(Signal(1:2000))
% ModemWave_down=downsample(ModemWave,6);
% figure;plot(ModemWave_down(1:2000))
% Signal_cut=Signal((pian_rece-Delay)*64:(pian_rece-Delay+15)*64);
% figure;plot(Signal_cut)
% figure;plot(Signal(1:1000))
%    figure;
%    plot(DemodWave(1:end))
%    DemodWave_down=downsample(DemodWave,2);
%       figure;
%    plot(DemodWave_down)
%          figure;
%    plot(DemodWave_down((pian_rece-5)*(cc/2)+20:(pian_rece-5+15)*(cc/2)+20))
%    figure;
%    plot(DemodWave((pian_rece-5)*64+40:(pian_rece-5+15)*64+40))
% %  RcvM_jiequ=RcvMatched((pian_rece)*16+10:(pian_rece+20)*16+10);
%    figure;
%    plot(RcvMatched((pian_rece)*64+40:(pian_rece+15)*64+40))
%       figure;
%    plot(RcvMatched)
%  figure;
%  plot(RcvMatched(10000:35000))
%  
% yI_r=[0 yI];
% yI_down=downsample(yI_r,cc/2);
% yI_jiequ=yI_down(pian_rece+2:pian_rece+2+20);
% figure;
% t=1:length(yI_down(pian_rece+2:pian_rece+2+15));
% scatter(t,yI_down(pian_rece+2:pian_rece+2+15))
% PP=yI_down(pian_rece+2:pian_rece+2+20);
% figure;
% plot(yI((pian_rece+2)*2:(pian_rece+2+20)*2))
% figure;plot(yI)
% figure;
% t=1:length(yI);
% scatter(t,yI)