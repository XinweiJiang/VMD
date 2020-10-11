% discussion on iteration convergence at different initial condition

clear;
clc;
close all;

rng(10001);             % random seed   
% rng('shuffle')
fs = 1024;          % sample frequency
t = (1:fs)/fs;    % discrete time step;
N = length(t);

fc = [30;100;300]; % frequency contents
s = sin(2*pi*fc*t);     % IMF

% B = randn(4,6);         % mixing matrix
% B = [-0.0664   -0.5757   -0.5445    0.4901
%     0.4575   -0.0874    0.2149    0.3464
%    -0.6186   -0.8098   -0.8087    0.6387
%    -0.6353    0.0715   -0.0580   -0.4815];
Br = randn(4,3);         % mixing matrix
x = Br*s;                % decomposition signals
x = awgn(x,20,'measured');

% solution parameters
alpha = 2000;
tau = 0;
K = 3;
DC = 0;
init = 4000;
tol = 1e-6;

figure
hold on
cout = 0;
for i=1:100
[u, u_hat, omega, Bc]=MVMD(x, alpha, tau, K, DC, init, tol);
if norm(fc-sort(omega(end,:)*fs)')<5
    for j=1:K
        if abs(omega(end,j)*fs-fc(1))<1
            plot(omega(:,j)*fs,0:size(omega,1)-1,'LineWidth',1.4, 'color','r')
        elseif abs(omega(end,j)*fs-fc(2))<1
            plot(omega(:,j)*fs,0:size(omega,1)-1,'LineWidth',1.4, 'color','g')
        else      
            plot(omega(:,j)*fs,0:size(omega,1)-1,'LineWidth',1.4, 'color','b') 
        end
    end
    cout = cout+1;
end
i
end


N = 15;
plot(ones(1,N)*30,0:N-1,'LineWidth',1.4, 'color','k','LineStyle',':')
plot(ones(1,N)*100,0:N-1,'LineWidth',1.4, 'color','k','LineStyle',':')
plot(ones(1,N)*300,0:N-1,'LineWidth',1.4, 'color','k','LineStyle',':')

xlabel('Frequency / Hz','fontsize',18)
ylabel('Iterations','fontsize',18)
set(gca,'FontSize',18);
xlim([0 420])
ylim([0 N])
box on
