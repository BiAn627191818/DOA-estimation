clc
clear all
close all


c=3e8;%光速
Tp=100e-6;%脉冲宽度
fp=1/Tp;%正交信号频率间隔
B=50e6;%带宽
u=B/Tp;%调频斜率
K=128;%码数
M=10;%发射阵元数
N=10;%接收阵元数
P=4;%目标个数
L=100;
f0=30e6;
w_length=c/f0;
dt=w_length/2;
dr=w_length/2;
DOD=pi*[-20,-15,-40,35]/180;%波离方向角
DOA=pi*[10,15,-20,20]/180;%波达方向角
fdp=[200,250,300,350];%多普勒频移


%DOD矩阵和DOA矩阵都是10*4,对于某个时刻

for p=1:P%快拍数
  Bt(:,p)=exp(j*2*pi*(0:M-1)*dt*sin(DOD(p))/w_length);
end

for p=1:P
    Ar(:,p)=exp(j*2*pi*(0:N-1)*dr*sin(DOA(p))/w_length);
end

%发射信号10*128

K=128;
HAR=hadamard(K);
S=HAR(1:10,:);

% RCS和多普勒频移组成对角阵4*4,与快拍l有关

X=zeros(N,K,L);
for l=1:L
    RCS=[2,5,9,20];
    Up=RCS.*exp(j*2*pi*fdp*Tp*l);
    Up=diag(Up);
    X(:,:,l)= Ar*Up*Bt.'*S+2*randn(N,K);
end

%匹配滤波
for l=1:L        
Yl(:,:,l)=X(:,:,l)*S'/K;
end

%匹配后滤波后，将矩阵变为M*N，L的矩阵

Y=zeros(M*N,L);
for l=1:L
    Y(:,l)=reshape(Yl(:,:,l).',M*N,1);
end

%求协方差阵与逆矩阵
Ry=1/L*Y*Y';
Ry_inv=inv(Ry);%

sita_t=(-90:0.5:90)*pi/180;
sita_r=(-90:0.5:90)*pi/180;
P_capon=zeros(length(sita_t),length(sita_r));
for m=1:length(sita_t)
    for n=1:length(sita_r)
        bt=exp(j*2*pi*dt*(0:M-1).'*sin(sita_t(m))/w_length);
        ar=exp(j*2*pi*dr*(0:N-1).'*sin(sita_r(n))/w_length);
        KLNK=kron(ar,bt);
        P_capon(m,n)=1/(KLNK'*Ry_inv*KLNK);
    end
end
P_capon=10*log10(abs(P_capon)/max(max(abs(P_capon))));
[X,Y]=meshgrid(-90:0.5:90);
    
mesh(X,Y,P_capon)

