% Copyright (c) 2020 Shuaishuai Liu. All rights reserved.

% We only permit to use these programs to verify our paper, "Multi-dimensional Variational Mode Decomposition and Its Short-time Counterpart".
% Other purposes are not permitted until further notice.

% complex signal decomposition

clear;
clc;
close all;

rng(1);             % random seed   
fs = 512;          % sample frequency
t = (1:2*fs)/fs;    % discrete time step;
N = length(t);

fc = [5;40;85;120;175]; % frequency contents
s = sin(2*pi*fc*t);     % IMF

Br = randn(4,5);         % real part of mixing matrix
Bi = randn(4,5);         % imaginary part of mixing matrix
B = Br+1j*Bi;
x = B*s;                % decomposition signals

% solution parameters
alpha = 2000;
tau = 0;
K = 5;
DC = 0;
init = 1;
tol = 1e-6;

% figure
% for i=1:100
x = [real(x);imag(x)];
[u, u_hat, omega, Bc]=MVMD(x, alpha, tau, K, DC, init, tol);
% plot(omega*fs,1:size(omega,1))
% hold on
% end
% hold off
% xlim([0 512])
% ylim([0 20])


figure
subplot(6,2,1)
plot(t,x(1:4,:))
ylabel('input')
title('Real(x)')
set(gca,'FontSize',12);

subplot(6,2,2)
plot(t,x(5:end,:))
title('Imag(x)')
set(gca,'FontSize',12);

for k=1:K
    subplot(6,1,k+1)
    if i<3
        plot(t,u(k,:),'LineWidth',1.4, 'color','k')
    else
        plot(t,u(k,:), 'color','k')
    end
    ylabel(['u_' int2str(k)])
    set(gca,'FontSize',12);
end
xlabel('Time / Sec')
set(gca,'FontSize',12);

% B
Bc = Bc(1:end/2,:)+1j*Bc((end/2+1):end,:);
for i = 1:K
    B(:,i) =   B(:,i)/B(1,i);
    Bc(:,i) =  Bc(:,i)/Bc(1,i);
end
% % (angle(B)-angle(Bc))./(angle(B)+eps*10)*100
% figure
% mac = MAC(B,Bc);
% bar3(mac)
% 
for i = 1:K
    B(:,i) =   B(:,i)/norm(B(:,i));
    Bc(:,i) =  Bc(:,i)/norm(Bc(:,i));
    acos(abs(B(:,i)'*Bc(:,i)))
end


