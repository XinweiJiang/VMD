% Copyright (c) 2020 Shuaishuai Liu. All rights reserved.

% We only permit to use these programs to verify our paper, "Multi-dimensional Variational Mode Decomposition and Its Short-time Counterpart".
% Other purposes are not permitted until further notice.

% long-term signal decomposition

clear;
clc;
close all;

% rng(7);             % random seed   
rng(995);             % random seed   
fs = 1024;          % sample frequency
t = (1:20*fs)/fs;    % discrete time step;
N = length(t);


s = [exp(0+0.1*cos(2*pi*0.15*t+pi/3)).*sin(2*pi*(100*t + 1.5*t.^2 + 2/pi*sin(0.9*pi*t)));
     exp(0+0.2*cos(2*pi*0.2*t+pi/4)).*sin(2*pi*(160*t - 1.5*t.^2 + 2/pi*sin(0.9*pi*t)));
     exp(0+0.2*cos(2*pi*0.1*t+pi/5)).*sin(2*pi*(190*t + 8/pi*sin(0.9*pi*t)));
     exp(0+0.3*cos(2*pi*0.2*t+pi/6)).*sin(2*pi*(200*t + 8/pi*sin(0.9*pi*t)))];     % IMF

fc = [100 + 3*t + 1.8*cos(0.9*pi*t);
      160 - 3*t + 1.8*cos(0.9*pi*t);
      190 + 7.2*cos(0.9*pi*t);
      200+ 7.2*cos(0.9*pi*t);]; 

% plot(t,fc)  

B = randn(3,4);         % mixing matrix
x = B*s;                % decomposition signals
% x = awgn(x,25,'measured');

% solution parameters
alpha = 2000;
tau = 1e-6*0;
K = 4;
DC = 0;
init = 1;
tol = 1e-6;
winLen = 512;
overlap = winLen-128;

[u, u_hat, omega, Bc]=STMVMD(x, alpha, tau, K, DC, init, tol, winLen, overlap);
omega = omega(:,[2 1 4 3]);
Bc = Bc(:,[2 1 4 3]);
u = u([2 1 4 3],:,:);

Nw = floor((length(x)-overlap)/(winLen-overlap));
index = ((0:Nw-1)*(winLen-overlap)+1)+((0:Nw-1)*(winLen-overlap)+winLen);
index = floor(index/2);
tc = t(index);

% figure
% plot(t,fc(1,:),'LineWidth',1.4, 'color','r')
% hold on
% plot(t,fc(2,:),'LineWidth',1.4, 'color','g')
% plot(t,fc(3,:),'LineWidth',1.4, 'color','b')
% plot(t,fc(4,:),'LineWidth',1.4, 'color','m')
% hold off
% xlabel('Time / Sec','fontsize',18)
% ylabel('Frequency / Hz','fontsize',18)
% set(gca,'FontSize',18);
% ylim([80 220])

figure
fcc = omega'*fs;
% plot(tc,fcc(1,:),'*','MarkerSize',3)
% hold on
% plot(tc,fcc(2,:),'*','MarkerSize',3)
% plot(tc,fcc(3,:),'*','MarkerSize',3)
% plot(tc,fcc(4,:),'*','MarkerSize',3)
plot(tc,fcc(1,:),'LineWidth',1.4, 'color','r')
hold on
plot(tc,fcc(2,:),'LineWidth',1.4, 'color','g')
plot(tc,fcc(3,:),'LineWidth',1.4, 'color','b')
plot(tc,fcc(4,:),'LineWidth',1.4, 'color','m')
hold off
xlabel('Time / Sec','fontsize',18)
ylabel('Frequency / Hz','fontsize',18)
set(gca,'FontSize',18);
ylim([80 220])


% stitching together

wc = 1 - (1:overlap)/overlap;           % linear weight of the current window
wn = 1 + ((-overlap):-1)/overlap;       % linear weight of the next window

% wc = 1./(1+exp(20*wn-10));            % nonlinear weight of the current window
% wn = 1./(1+exp(-20*wn+10));           % nonlinear weight of the next window

wc = repmat(wc,K,1);
wn = repmat(wn,K,1);

Ns = (Nw-1)*(winLen-overlap)+winLen;
T = t(1:Ns);
U = zeros(K,Ns);

U(:,1:winLen) = squeeze(u(:,:,1));
for nw = 1:Nw-1
    Uc = squeeze(u(:,:,nw));
    Un = squeeze(u(:,:,nw+1));
    Unc = [Uc(:,1:winLen-overlap) (wc.*Uc(:,winLen-overlap+1:end) ...
            +wn.*Un(:,1:overlap)) Un(:,overlap+1:end)];
    U(:,((nw-1)*(winLen-overlap)+floor(winLen/2)):((nw-1)*(winLen-overlap)+winLen)) = ...
         Unc(:,floor(winLen/2):winLen);
end
U(:,end-floor(winLen/2)+1:end) = Un(:,end-floor(winLen/2)+1:end);

% U = U([2 1 4 3],:);
figure
subplot(5,1,1)
plot(t,x)
fontsize = 12;
ylabel('input','fontsize',fontsize)
set(gca,'FontSize',fontsize);
for k=1:K
    subplot(5,1,k+1)
    plot(T,U(k,:),'color','k')
    ylabel(['u_' int2str(k)],'fontsize',fontsize)
    ylim([-3 3])
    set(gca,'FontSize',fontsize);
%     xlim([8 12])
end
xlabel('Time / Sec','fontsize',fontsize)
set(gca,'FontSize',fontsize);


% Uc = squeeze(u(:,:,1));
% Un = squeeze(u(:,:,2));
% Unc = [Uc(:,1:winLen-overlap) (wc.*Uc(:,winLen-overlap+1:end) +wn.*Un(:,1:overlap)) Un(:,overlap+1:end)];
% figure
% ch = 1;
% subplot(2,2,1)
% % plot(T(1:640),[Uc(ch,:) zeros(1,128)])
% plot(T(1:512),Uc(ch,:))
% xlim([0.1 0.15])
% subplot(2,2,2)
% % plot(T(1:640),[zeros(1,128) Un(ch,:)])
% plot(T(129:640),Un(ch,:))
% xlim([0.1 0.15])
% subplot(2,1,2)
% plot(T(1:640),U(ch,1:640))
% xlim([0.1 0.15])


% axes('Position',[0.45,0.45,0.1,0.09]); % Éú³É×ÓÍ¼   
% plot(tc,fcc(1,:),'LineWidth',1.4, 'color','r')
% hold on
% plot(tc,fcc(2,:),'LineWidth',1.4, 'color','g')   
% xlim([9 11])
% ylim([120 140])


% figure
% [S, F, T] = spectrogram(x(3,:),hanning(winLen),overlap,winLen,fs,'oneside');
% imagesc(T, F, abs(S));
% set(gca,'YDir','normal'); 
% xlabel('Time / Sec','fontsize',18)
% ylabel('Frequency / Hz','fontsize',18)
% set(gca,'FontSize',18);
% % xticks([0 5 10 15 20])
% ylim([80 220])
