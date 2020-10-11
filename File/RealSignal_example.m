% Copyright (c) 2020 Shuaishuai Liu. All rights reserved.

% We only permit to use these programs to verify our paper, "Multi-dimensional Variational Mode Decomposition and Its Short-time Counterpart".
% Other purposes are not permitted until further notice.

% real signal decomposition

clear;
clc;
close all;

rng(1);             % random seed   
fs = 1024;          % sample frequency
t = (1:fs)/fs;    % discrete time step;
N = length(t);

fc = [2;24;36]; % frequency contents
s = sin(2*pi*fc*t);     % IMF


Br = randn(2,3);         % mixing matrix
% Br(1,2) = 0;
% Br(2,1) = 0;
x = Br*s;                % decomposition signals
% x = awgn(x,20,'measured');

% solution parameters
alpha = 2000;
tau = 0;
K = 3;
DC = 0;
init = 1;
tol = 1e-6;

tic
[u, u_hat, omega, Bc]=MVMD(x, alpha, tau, K, DC, init, tol);
toc

figure
subplot(4,2,1)
plot(t,x(1,:),'LineWidth',1.2 ,'color','k')
ylim([-3 3])
ylabel('input')
title('Channel x_1')
set(gca,'FontSize',12);

subplot(4,2,2)
plot(t,x(2,:),'LineWidth',1.2 ,'color','k')
ylim([-3 3])
title('Channel x_2')
set(gca,'FontSize',12);

for i=1:3
    subplot(4,2,2*i+1)
    plot(t,Bc(1,i)*u(i,:),'LineWidth',1.4, 'color','k')
    ylabel(['u_' int2str(i)])
    if i<=2
        ylim([-2 2])
    end
    set(gca,'FontSize',12);
end
xlabel('Time / Sec')
set(gca,'FontSize',12);

for i=1:3
    subplot(4,2,2*i+2)
    plot(t,Bc(2,i)*u(i,:),'LineWidth',1.4, 'color','k')
    if i<=2
        ylim([-2 2])
    end
    set(gca,'FontSize',12);
end
xlabel('Time / Sec')
set(gca,'FontSize',12);

% cosbi = zeros(1,3);
% for i = 1:K
%     Br(:,i) =   Br(:,i)/norm(Br(:,i));
%     Bc(:,i) =  Bc(:,i)/norm(Bc(:,i));
%     cosbi(i) = Br(:,i)'*Bc(:,i);
% end
% cosbi

% figure
% mac = MAC(Br,Bc);
% bar3(mac)
% 
% figure
% for k=1:K
%     subplot(2,2,k)
%     plot(t,u(k,:))
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
