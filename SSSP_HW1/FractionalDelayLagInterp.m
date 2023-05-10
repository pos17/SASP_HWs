%------------------------------------------%
%        *** SSSP - HOMEWORK #3 ***        %
%------------------------------------------%
%     Fractional delay using lagrange      %
%               Interpolation              %
%                                          % 
%------------------------------------------%
% Name:                                    %
% Student ID:                              %
%------------------------------------------%
clear
close all
clc

D = 3.3; %tau*fs
fs = 44100;

%% tests with N = 1 and N = 2 
% implemented to test existence conditions
%
% N = 1;
% 
% h_int_1 = LagrangeInterpolation(N,D);
% n_1 = 0:N;
% 
% 
% [H_int_1, f_1] = freqz(h_int_1, 1, 'half', fs);
% 
% 
% N = 2;
% 
% h_int_2 = LagrangeInterpolation(N,D);
% n_2 = 0:N;
% [H_int_2, f_2] = freqz(h_int_2, 1, 'half', fs);
% 
N = 4;

h_int_4 = LagrangeInterpolation(N,D);
n_4= 0:N;
[H_int_4, f_4] = freqz(h_int_4, 1, 'half', fs);

N = 5;

h_int_5 = LagrangeInterpolation(N,D);
n_5= 0:N;
[H_int_5, f_5] = freqz(h_int_5, 1, 'half', fs);

% optimal filter settings
N = 2*floor(D);
n_opt = 0:N;
h_int_opt = LagrangeInterpolation(N,D);
[H_int_opt, f_opt] = freqz(h_int_opt, 1, 'half', fs);

% optimal filter settings
N = 7
n_7 = 0:N;
h_int_7 = LagrangeInterpolation(N,D);
[H_int_7, f_7] = freqz(h_int_7, 1, 'half', fs);

% optimal filter settings
N = 8%4*floor(D);
n_8 = 0:N;
h_int_8 = LagrangeInterpolation(N,D);
[H_int_8, f_8] = freqz(h_int_8, 1, 'half', fs);

% optimal filter settings
N = 9%5*floor(D);
n_9 = 0:N;
h_int_9 = LagrangeInterpolation(N,D);
[H_int_9, f_9] = freqz(h_int_9, 1, 'half', fs);



%% Lagrange interpolation plots


figure('Renderer', 'painters', 'Position', [10 10 1000 800])
% stem(n_1, h_int_1),
% hold on 
% stem(n_2, h_int_2),
% hold on 
stem(n_4, h_int_4),
hold on 
stem(n_5, h_int_5),
hold on 
stem(n_opt, h_int_opt),
hold on 
stem(n_7, h_int_7),
hold on 
stem(n_8, h_int_8),
hold on 
stem(n_9, h_int_9),
xlabel('Samples',Interpreter="latex"),
ylabel('Impulse response',Interpreter="latex")
title('Lagrange interpolation',Interpreter="latex")
%legend(["N=1","N=2","N=3","N= floor(2*D)","N=7","N= floor(4*D)","N= floor(5*D)"],'Location','northwest',Interpreter="latex");
legend(["N=4","N=5","N=6","N=7","N=8","N=9"],'Location','northwest',Interpreter="latex");
saveas(gcf,"Plots/FiltersIR.png")
%[H_int, f] = freqz(h_int_1, 1, 'half', fs);

figure('Renderer', 'painters', 'Position', [10 10 1000 1200])

subplot(3,1,1)
subtitle('Frequency response',Interpreter="latex");
% plot(f_1*fs/pi/2, abs(H_int_1));
% hold on
% plot(f_2*fs/pi/2, abs(H_int_2));
% hold on
%norm_freq = (0,1,f_4*fs/pi/2);
plot(f_4/(pi), abs(H_int_4));
hold on
plot(f_5/(pi), abs(H_int_5));
hold on
plot(f_opt/(pi), abs(H_int_opt));
hold on
plot(f_7/(pi), abs(H_int_7));
hold on
plot(f_8/(pi), abs(H_int_8));
hold on
plot(f_9/(pi), abs(H_int_9));
ylim([0 1.7]);
%xlim([0 fs/2]);
xlabel('Frequency [Hz]',Interpreter="latex"), ylabel('Magnitude [dB]',Interpreter="latex");
%legend(["N=1","N=2","N=3","N= floor(2*D)","N=7","N= floor(4*D)","N= floor(5*D)"],'Location','northwest',Interpreter="latex");
legend(["N=4","N=5","N=6","N=7","N=8","N=9"],'Location','northwest',Interpreter="latex");
hold on
%axis([0, fs/2, 0 2]),
grid

subplot(3,1,2)
subtitle('Group Delay',Interpreter="latex");
% grpdelay(h_int_1, 1, 1024, 'half', fs);
% hold on 
% grpdelay(h_int_2, 1, 1024, 'half', fs);
% hold on 
grpdelay(h_int_4, 1, 1024);
hold on
grpdelay(h_int_5, 1, 1024);
hold on 
grpdelay(h_int_opt, 1, 1024);
hold on 
grpdelay(h_int_7, 1, 1024);
hold on 
grpdelay(h_int_8, 1, 1024);
hold on 
grpdelay(h_int_9, 1, 1024);
ylim([1 4]);
xlabel('Normalized Frequency $(x \pi [rad/sample])$',Interpreter="latex"), ylabel('Group Delay [samples]',Interpreter="latex");
%xlim([0 fs/2]);
%legend(["N=1","N=2","N=3","N= floor(2*D)","N=7","N= floor(4*D)","N= floor(5*D)"],'Location','northwest',Interpreter="latex");
legend(["N=4","N=5","N=6","N=7","N=8","N=9"],'Location','northwest',Interpreter="latex");

subplot(3,1,3)
subtitle('Group Delay',Interpreter="latex");
phasedelay(h_int_4, 1, 1024);
hold on 
phasedelay(h_int_5, 1, 1024);
hold on 
phasedelay(h_int_opt, 1, 1024);
hold on 
phasedelay(h_int_7, 1, 1024);
hold on 
phasedelay(h_int_8, 1, 1024);
hold on 
phasedelay(h_int_9, 1, 1024);
xlabel('Normalized Frequency $(x \pi [rad/sample])$',Interpreter="latex"), ylabel('Phase Delay [samples]',Interpreter="latex");

%xlim([0 1024]);
%legend(["N=1","N=2","N=3","N= floor(2*D)","N=7","N= floor(4*D)","N= floor(5*D)"],'Location','northwest',Interpreter="latex");
legend(["N=4","N=5","N=6","N=7","N=8","N=9"],'Location','northwest',Interpreter="latex");
sgtitle('Filter response orders comparison',Interpreter="latex");
saveas(gcf,"Plots/FreqResponseComparison.png")


x = linspace(0,2*floor(D),100);
sinc_y = sinc(x-(D));

figure('Renderer', 'painters', 'Position', [10 10 1000 800])
stem(n_opt, h_int_opt),
hold on 
plot(x,sinc_y);