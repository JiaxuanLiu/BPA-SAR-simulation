%%%%%%%%%%%%%%%%%% 声明 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% 该代码使用BP成像方式，基于LFM波体制 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% 该代码仅限个人用途，以研究为主 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% 改版增加了注释 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% 修改参数时需要仔细考虑，这一版为远场成像 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% 模拟的是飞机雷达对地面目标成像 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% 2021.7.8 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Made by JiaxuanLiu %%%%%%%%%%%%%%%%%%

%%  参数设置
clc,clear all,close all
C = 3e8;                                                % 光速
fc = 5.3e9;                                             % 中心频率
lambda = C/fc;                                          % 波长
D = 4;                                                  % 孔径大小
V = 150;                                                % 雷达移动速度
Kr = 20e12;                                             % 调频斜率
Tr = 2.5e-6;                                            % 一个chirp升频率时间
sq_ang = 3.5/180*pi;                                    % 雷达探测区域角度

Br = Kr*Tr;                                             % 信号带宽
Frfactor = 1.2;                                         % 采样倍数（距离向）
Fr = Br*Frfactor;                                       % 采样频率（距离向）
Ba = 0.886*2*V*cos(sq_ang)/D;                           % 方位向信号带宽
Fafactor = 1.2;                                         % 采样倍数（方位向）
Fa = Ba*Fafactor;                                       % 采样频率（方位向）

R_near = 2e4;                                           % 成像区域近点
R_far = R_near+1000;                                    % 成像区域远点

% 以下为仿真场景设置
La_near = 0.886*R_near*lambda/cos(sq_ang)^2/D;
La_far = 0.886*R_far*lambda/cos(sq_ang)^2/D;
Tc_near = -R_near*tan(sq_ang)/V;
Tc_far = -R_far*tan(sq_ang)/V;
fdc = 2*V*sin(sq_ang)/lambda;
Y_min = V*Tc_far;
Y_max = Y_min+100;

Rmin = sqrt(R_near^2+(Tc_near*V+La_near/2)^2);
Rmax = sqrt(R_far^2+(Tc_far*V-La_far/2)^2);

%% 回波设置
Nr = (2*Rmax/C-2*Rmin/C+Tr)*Fr;
Nr = 2^nextpow2(Nr);
tr = linspace(-Tr/2+2*Rmin/C,Tr/2+2*Rmax/C,Nr);         % 快时间
Fr = (Nr-1)/(Tr/2+2*Rmax/C-(-Tr/2+2*Rmin/C));
Na = ((Tc_near+La_near/2/V)-(Tc_far-La_far/2/V))*Fa;
Na = 2^nextpow2(Na);
ta = linspace(Tc_far-La_far/2/V,Tc_near+La_near/2/V,Na);% 满时间
Fa = (Na-1)/(Tc_near+La_near/2/V-(Tc_far-La_far/2/V));

% 方针目标设置
Rpt = [R_near R_near+500 R_near+1000];
Ypt = [0 0 0];
La = 0.886*Rpt*lambda/cos(sq_ang)^2/D;
Tc = -Rpt*tan(sq_ang)/V;
Npt = length(Rpt);

Y_high = max(Ypt)+50;
Y_low = min(Ypt)-50;
% R_left = R_near-50;
% R_right = R_far+50;
R_left = Rmin;
R_right = Rmax;
row    = tr*C/2;                        %列
col    = ta*V;    
% 回波仿真
sig = zeros(Na,Nr);
for k = 1:Npt
    delay = 2/C*sqrt(Rpt(k)^2+(Ypt(k)-ta*V).^2);
    Dr = ones(Na,1)*tr-delay'*ones(1,Nr);
    sig = sig+exp(1j*pi*Kr*Dr.^2-1j*2*pi*fc*delay'*ones(1,Nr))...
        .*(abs((ta-Ypt(k)/V-Tc(k))'*ones(1,Nr))<=La(k)/2/V).*(abs(Dr)<=Tr/2);
end
figure,imagesc(row,col,real(sig)) ;caxis([0 1]);axis xy

%% 距离压缩
sig_rd = fft(sig,[],2);
fr = -1/2:1/Nr:(1/2-1/Nr);
fr = fftshift(fr*Fr);
filter_r = ones(Na,1)*exp(1j*pi*fr.^2/Kr);
sig_rd = sig_rd.*filter_r;
nup = 10;
Nr_up = Nr*nup;
nz = Nr_up-Nr;
dtr = 1/nup/Fr;
sig_rd_up = zeros(Na,Nr_up);
sig_rd_up = [sig_rd(:,1:Nr/2),zeros(Na,nz),sig_rd(:,(Nr/2+1):Nr)];
sig_rdt = ifft(sig_rd_up,[],2);
figure,imagesc(abs(sig_rdt));

%% 网格化
R = zeros(1,Nr);
for ii = 1:Nr
    R(1,ii) = R_left+(R_right-R_left)/(Nr-1)*(ii+1);
end
Y = zeros(1,Na);
for ii = 1:Na
    Y(1,ii) = Y_low+(Y_high-Y_low)/(Na-1)*(ii-1);
end
R = ones(Na,1)*R;
Y = Y'*ones(1,Nr);

%% 成像
f_back = zeros(Na,Nr);
h = waitbar(0,'BPA');
for ii = 1:Na
    R_ij = sqrt(R.^2+(Y-V*ta(ii)).^2);
    t_ij = 2*R_ij/C;
    t_ij = round((t_ij-(2*Rmin/C-Tr/2))/dtr);       %时域转换成时域点数
    it_ij = (t_ij>0&t_ij<=Nr_up);
    t_ij = t_ij.*it_ij+Nr_up*(1-it_ij);
    sig_rdta = sig_rdt(ii,:);
    sig_rdta(Nr_up) = 0;
    f_back = f_back+sig_rdta(t_ij).*exp(1j*4*pi*R_ij/lambda);       %（球面）累加曲线
    waitbar(ii/Na);
    ii
end
close(h);
figure,imagesc(abs(f_back)),suptitle('升采样率为10时成像结果'),xlabel('距离维'),ylabel('方位维');