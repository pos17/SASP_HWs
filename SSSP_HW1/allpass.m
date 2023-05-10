%% First order shelving filter
%  Lucio Bianchi - 03/04/2015

clear
close all
clc

%% Parameters
fs = 44100; % sampling frequency [Hz]
Nfft = 4096; % number of FFT bins
f = 0:fs/Nfft:fs-fs/Nfft;

%% Low frequency shelving
fc = 2000; % cutoff frequency [Hz]
Wc = 2*pi*fc/fs; % normalized cutoff frequency

% Boost
G = 20; % logarithmic gain [dB]
V0 = 10^(G/20);
H0 = V0-1;

a1 = -0.95;

num = [a1^4,4*a1^3,6*a1^2,4*a1,1];
den = [1,4*a1, 6*a1^2,4*a1^3,a1^4];
zplane(num,den);
H = freqz(num, den, f(1:Nfft/2+1), fs);
saveas(gcf,'Plots/zplane',"png");

figure,
subplot(2,1,1)
semilogx(f(1:Nfft/2+1), 20*log10(abs(H))), grid, title(['Low frequency boost, G = ' num2str(G) ' dB'])
xlabel('Frequency [Hz]'), ylabel('Magnitude [dB]')
axis([20 f(Nfft/2+1) -20 20])

subplot(2,1,2)
plot(f(1:Nfft/2+1), angle(H));
ylabel('\angle H'); % label left y-axis
xlabel('\omega'); % label x-axis
title('TF - Phase');
grid on

