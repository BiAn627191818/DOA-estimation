clc
clear all
close all


c=3e8;%����
Tp=100e-6;%������
fp=1/Tp;%�����ź�Ƶ�ʼ��
B=50e6;%����
u=B/Tp;%��Ƶб��
K=128;%����
M=10;%������Ԫ��
N=10;%������Ԫ��
P=4;%Ŀ�����
L=100;
f0=30e6;
w_length=c/f0;
dt=w_length/2;
dr=w_length/2;
DOD=pi*[-20,-15,-40,35]/180;%���뷽���
DOA=pi*[10,15,-20,20]/180;%���﷽���
fdp=[200,250,300,350];%������Ƶ��


%DOD�����DOA������10*4,����ĳ��ʱ��

for p=1:P%������
  Bt(:,p)=exp(j*2*pi*(0:M-1)*dt*sin(DOD(p))/w_length);
end

for p=1:P
    Ar(:,p)=exp(j*2*pi*(0:N-1)*dr*sin(DOA(p))/w_length);
end

%�����ź�10*128

K=128;
HAR=hadamard(K);
S=HAR(1:10,:);

% RCS�Ͷ�����Ƶ����ɶԽ���4*4,�����l�й�

X=zeros(N,K,L);
for l=1:L
    RCS=[2,5,9,20];
    Up=RCS.*exp(j*2*pi*fdp*Tp*l);
    Up=diag(Up);
    X(:,:,l)= Ar*Up*Bt.'*S+2*randn(N,K);
end

%ƥ���˲�
for l=1:L        
Yl(:,:,l)=X(:,:,l)*S'/K;
end

%ƥ����˲��󣬽������ΪM*N��L�ľ���

Y=zeros(M*N,L);
for l=1:L
    Y(:,l)=reshape(Yl(:,:,l).',M*N,1);
end

%��Э�������������
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

