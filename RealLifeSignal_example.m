% real-life signal decomposition

clear;
clc;
close all;

% 数据源：https://dataverse.tdl.org/dataset.xhtml?persistentId=doi:10.18738/T8/SS2NHB
% EEG_Cat_Study4_II_II_S1.bdf
% 提取时间：10~30s
% 提取通道：35 36 37通道AF8/AF4/AFz
data = load('eeg.mat');
x = data.y(2:end,:);
fs = 256;
t = 10:1/fs:(20-1/fs);
Ns = max(size(t));
x = x(1:Ns,:);

% solution parameters
alpha = 2000;
tau = 0;
K = 4;
DC = 0;
init = 1;
tol = 1e-6;

tic
[u, u_hat, omega, Bc]=MVMD(x, alpha, tau, K, DC, init, tol);
toc

FontSize = 12;
figure
subplot(5,1,1)
hold on
plot(t,x(:,1),'LineWidth',1.4, 'color','r')
plot(t,x(:,2),'LineWidth',1.4, 'color','g','LineStyle','--')
plot(t,x(:,3),'LineWidth',1.4, 'color','b','LineStyle','-.')
hold off
box on
ylabel('input')
set(gca,'FontSize',FontSize);

u = u([1 2 4 3],:);

for k=1:K
    subplot(5,1,k+1)
    plot(t,u(k,:),'LineWidth',1.4, 'color','k')
    ylabel(['u_' int2str(k)])
    if k==1
%        ylim([-100 200]) 
    end
    set(gca,'FontSize',FontSize);
end
xlabel('Time / Sec')
set(gca,'FontSize',FontSize);



figure
% subplot(2,2,1)
% semilogy((0:Ns-1)/Ns*fs,abs(fft(x)/Ns).^2)
% xlim([0 70])
% ylim([1e-5 1e2])
for k=1:K
    subplot(2,2,k)
    [Pxx,F] = pwelch(u(k,:),[],[],[],fs,'onesided');
   
    plot(F,10*log10(Pxx),'LineWidth',1.4,'color','k')
%     ylabel(['u_' int2str(k)])
    if k==1 || k==3
        ylabel(['Power (dB)'])
    end
    xlim([0 70])
%     ylim([-inf 1e2])
    xlabel('Freqency / Hz')
    set(gca,'FontSize',FontSize);  
end

set(gca,'FontSize',FontSize);

% plot(x)
% 
% winLen = 256;
% overlap = winLen-128;
% figure
% spectrogram(x(:,3),hanning(winLen),overlap,winLen,fs,'oneside');
% [S, F, T] = spectrogram(x(:,3),hanning(winLen),overlap,winLen,fs,'oneside');
% imagesc(T, F, abs(S));
% set(gca,'YDir','normal'); 
% xlabel('Time / Sec','fontsize',18)
% ylabel('Frequency / Hz','fontsize',18)
% set(gca,'FontSize',18);
% % xticks([0 5 10 15 20])
% ylim([80 220])

