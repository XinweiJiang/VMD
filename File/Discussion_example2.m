% Copyright (c) 2020 Shuaishuai Liu. All rights reserved.

% We only permit to use these programs to verify our paper, "Multi-dimensional Variational Mode Decomposition and Its Short-time Counterpart".
% Other purposes are not permitted until further notice.

% discussion on CPU consume

clear;
clc;
close all;

rng(1);             % random seed   
fs = 1024;          % sample frequency
t = (1:fs)/fs;    % discrete time step;
N = length(t);

fc = [30;100;200;300;400]; % frequency contents
s = sin(2*pi*fc*t);     % IMF


% solution parameters
alpha = 2000;
tau = 0;
K = 5;
DC = 0;
init = 1;
tol = 1e-6;

t = zeros(3,63);

for i = 2:64
i    
Br = randn(i,5);         % mixing matrix
x = Br*s;                % decomposition signals
% x = awgn(x,15,'measured');

tic
for ch=1:size(x,1)
    [u, u_hat, omega] = VMD(x(ch,:), alpha, tau, K, DC, init, tol);
end
% toc 
t(1,i-1) = toc;

tic
for ch=1:size(x,1)
    [u, u_hat, omega] = MVVMD(x, alpha, tau, K, DC, init, tol);
end
% toc
t(2,i-1) = toc;

tic
[u, u_hat, omega, Bc]=MVMD(x, alpha, tau, K, DC, init, tol);
% toc
t(3,i-1) = toc;

end

semilogy(2:64,t(1,:),'color','r','LineWidth',1.4,'lineStyle','-')
hold on
semilogy(2:64,t(2,:),'color','b','LineWidth',1.4,'lineStyle','-')
semilogy(2:64,t(3,:),'color','m','LineWidth',1.4,'lineStyle','-')
hold off
% hl = legend('VMD','MVVMD','MDVMD');
% set(hl,'box','off')

xlabel('Dimension','fontsize',12)
ylabel('Computional time / Sec','fontsize',12)

set(gca,'FontSize',18);
xlim([2 64])

% save('CPUtime_SNR00dB.mat','t')


% figure
% mac = MAC(Br,Bc);
% bar3(mac)

% figure
% for k=1:K
%     subplot(2,2,k)
%     plot(u(k,:))
% end
% 
% figure
% for k=1:K
%     subplot(2,2,k)
%     plot((0:N-1)/N*fs,abs(fft(u(k,:))))
%     xlim([0 fs/2])
% end
% 
% figure
% plot(omega*fs,1:size(omega,1))
