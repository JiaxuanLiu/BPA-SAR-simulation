clear all;close all;clc;

Fc = 25.5e9;
C = physconst('lightspeed');
Band = 4e9;
TimeUp = 640e-6;
TimeTotal = 1000e-6;
Fs = 25e5;
SamplePerChripUp = double(int32(TimeUp*Fs));
SamplePerChrip = double(int32(TimeTotal*Fs));
Nr = SamplePerChripUp;
Tfast = (0:SamplePerChripUp-1)/SamplePerChripUp*TimeUp;
Kr = Band/TimeUp;
dr = C/2/Band;
lambda = C/Fc;
theta = 120;
RadarVelocity = 5;
Timer = 4;

Na = Timer/TimeTotal;
RadarPos = [(-Timer/2:TimeTotal:Timer/2-1/Fs)*RadarVelocity;zeros(1,Na);zeros(1,Na)];
TargetPos = [0,10,0];
Range = sqrt((RadarPos(1,:)-TargetPos(1)).^2+(RadarPos(2,:)-TargetPos(2)).^2+(RadarPos(3,:)-TargetPos(3)).^2);
isFull = acosd(TargetPos(2)./Range)<theta/2;
% figure,plot(isFull);
N_up = 10;
tau = 2*Range/C;
Srnm = exp(1j*2*pi*(Fc*(Tfast-tau.')+1/2*Kr*(Tfast-tau.').^2)).*isFull.';
SrnmOri = exp(1j*2*pi*(Fc*Tfast+1/2*Kr*Tfast.^2));
SrnmFs = conj(Srnm).*SrnmOri;
rwin = hamming(SamplePerChripUp)*ones(1,Na);
SrnmFs = SrnmFs.*rwin.';
SrnmFt = fft(SrnmFs,SamplePerChripUp*N_up,2);
dr = dr/N_up;
Nr = Nr*N_up;
figure,plot((0:SamplePerChripUp*N_up-1)*dr,abs(SrnmFt(1000,:)));

%%
L = 0.1;
pix = 0.01;
x = -50*L:pix:50*L;
y = (-50*L:pix:50*L)+TargetPos(2);
[X,Y] = meshgrid(x,y);
img = zeros(size(x,2),size(y,2));
h = waitbar(0,'正在BP计算');
for ix = 1:Na
    RadarPosNow = RadarPos(:,ix);
    Ran = sqrt((RadarPosNow(1)-X).^2+(RadarPosNow(2)-Y).^2+RadarPosNow(3).^2);
    Ran = Ran.*(acosd(Y./Ran)<=60);
    deltaN = round(Ran/dr) + 1;
    deltaN_ij = (deltaN > 0 & deltaN <= Nr);
    deltaN = deltaN.*deltaN_ij + Nr*(1 - deltaN_ij);
    sig_rdta=SrnmFt(ix,:);
    sig_rdta(Nr)=0;
    img=img+sig_rdta(deltaN).*exp(1i*4*pi/lambda.*Ran);
    waitbar(ix/Na);
end
close(h);
figure,imagesc(abs(img));
