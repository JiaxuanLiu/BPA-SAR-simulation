%%%%%%%%%%%%%%%%%% 声明 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% 该代码使用BP成像方式，基于FMCW波体制 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% 该代码仅限个人用途，以研究为主 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% 改版增加了注释 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% 修改参数时需要仔细考虑，我这里算是近场成像 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% 所以会出现六芒星 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% 2021.7.7 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Made by JiaxuanLiu %%%%%%%%%%%%%%%%%%


clear all;close all;clc;

% 参数
Fc = 25.5e9;                                                    % 中心频率
C = physconst('lightspeed');                                    % 光速
Band = 4e9;                                                     % 带宽
TimeUp = 640e-6;                                                % 升频率时间
TimeTotal = 1000e-6;                                            % 一个chirp周期总时间
Fs = 25e5;                                                      % 采样频率
SamplePerChripUp = double(int32(TimeUp*Fs));                    % 一个chirp周期升频率采样点数
SamplePerChrip = double(int32(TimeTotal*Fs));                   % 一个chirp周期总采样点数
Nr = SamplePerChripUp;                                          % Nr
Tfast = (0:SamplePerChripUp-1)/SamplePerChripUp*TimeUp;         % 快时间序列
Kr = Band/TimeUp;                                               % 调频率
dr = C/2/Band;                                                  % 距离分辨率
lambda = C/Fc;                                                  % 信号波长
theta = 120;                                                    % 雷达有效探照角度
RadarVelocity = 5;                                              % 雷达运动速度
Timer = 4;                                                      % 仿真时间

Na = Timer/TimeTotal;                                           % 方位向位置总点数
RadarPos = [(-Timer/2:TimeTotal:Timer/2-1/Fs)*RadarVelocity; ...% 雷达位置
    zeros(1,Na);zeros(1,Na)];
TargetPos = [0,10,0];                                           % 仿真点位置
Range = sqrt((RadarPos(1,:)-TargetPos(1)).^2+(RadarPos(2,:)-... % 雷达与仿真点距离
    TargetPos(2)).^2+(RadarPos(3,:)-TargetPos(3)).^2);
isFull = acosd(TargetPos(2)./Range)<theta/2;                    % 判断目标是否在探照面内

N_up = 10;                                                      % 升采样倍数
tau = 2*Range/C;                                                % 延迟时间
Srnm = exp(1j*2*pi*(Fc*(Tfast-tau.')+1/2*Kr* ...                % 回波
    (Tfast-tau.').^2)).*isFull.';
SrnmOri = exp(1j*2*pi*(Fc*Tfast+1/2*Kr*Tfast.^2));              % 参考信号
SrnmFs = conj(Srnm).*SrnmOri;                                   % 得到中频信号
rwin = hamming(SamplePerChripUp)*ones(1,Na);                    % 汉明窗
SrnmFs = SrnmFs.*rwin.';                                        % 加窗降噪
SrnmFt = fft(SrnmFs,SamplePerChripUp*N_up,2);                   % 距离维FFT
dr = dr/N_up;                                                   % 更新距离分辨率
Nr = Nr*N_up;                                                   % 更新Nr
figure,plot((0:SamplePerChripUp*N_up-1)*dr,abs(SrnmFt(1000,:)));% 随便找个点看看

%% 
% 开始构建成像区域
L = 0.1;                                                        % 基本距离单元（随便设得）
pix = 0.01;                                                     % 像素距离
x = -50*L:pix:50*L;                                             % 构建x轴，我这里是方位维
y = (-50*L:pix:50*L)+TargetPos(2);                              % 以目标为中心构建y轴，对应距离维
[X,Y] = meshgrid(x,y);                                          % 构建
img = zeros(size(x,2),size(y,2));                               % 预构图像
h = waitbar(0,'正在BP计算');                                    % 加个进度条
for ix = 1:Na
    RadarPosNow = RadarPos(:,ix);                               % 找出现在的雷达位置
    Ran = sqrt((RadarPosNow(1)-X).^2+(RadarPosNow(2)-Y).^2 ...  % 计算当前时刻雷达和成像区域每一个点的距离
        +RadarPosNow(3).^2);
    Ran = Ran.*(acosd(Y./Ran)<=60);                             % 只考虑雷达照射范围内的点
    deltaN = round(Ran/dr) + 1;                                 % 计算距离对应的点数
    deltaN_ij = (deltaN > 0 & deltaN <= Nr);                    % 考虑合理区间
    deltaN = deltaN.*deltaN_ij + Nr*(1 - deltaN_ij);            % 考虑合理区间
    sig_rdta=SrnmFt(ix,:);                                      % 找出当前时刻雷达采集到的数据
    sig_rdta(Nr)=0;                                             % 找出合理区间
    img=img+sig_rdta(deltaN).*exp(1i*4*pi/lambda.*Ran);         % 将成像区域每个点在该时刻雷达数据中的点进行对应
    waitbar(ix/Na);                                             % 更新进度条
end
close(h);                                                       % 关闭进度条
figure,imagesc(abs(img));                                       % 成像看看
