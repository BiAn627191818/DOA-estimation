%%%%
%%均匀线阵下PM的DOA估计

clear;
M = 8;
lambda = 1;
d = 0.5;
L = 200;
theta = [-35 0 30 45];
K = length(theta);
snr = 5;
rad = pi/180;
A = exp(-j*2*pi*[0:M-1].'*d*sin(theta*rad)/lambda);
s = 10^(snr/10)*sign(2*(rand(K,L))-1);
N = 1/sqrt(2)*(randn(M,L)+j*randn(M,L));
X = A*s+N;
Rx = X*X'/L;
G = Rx(:,1:K);
H = Rx(:,K+1:end);
P = inv(G'*G)*G'*H;
Q = [P;-eye(M-K)];
afa = -90:0.1:90;
sp = zeros(1,length(afa));
for i = 1:length(afa)
    a = exp(-j*2*pi*[0:M-1].'*d*sin(afa(i)*rad)/lambda);
    sp(i) = 1/(a'*(Q*Q')*a);
end
sp = 10*log10(abs(sp));
plot(afa,sp);