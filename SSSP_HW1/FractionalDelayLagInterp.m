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

zplane(h_int_5,1)

%% Lagrange interpolation plots

lw = 1
fontSize = 15
figure('Renderer', 'painters', 'Position', [10 10 1000 1500])

subplot(3,2,1)
stem(n_4, h_int_4,LineWidth=lw),
hold on 
x = linspace(0,4,100);
sinc_y = sinc(x-(D));
plot(x,sinc_y,LineWidth=lw);
xlabel('Samples',Interpreter="latex",FontSize=fontSize),
ylabel('Impulse response',Interpreter="latex",FontSize=fontSize)
xlim([-1 5])
grid
title('N = 4',Interpreter="latex",FontSize=fontSize)


subplot(3,2,2)
stem(n_5, h_int_5,LineWidth=lw),
x = linspace(0,5,100);
sinc_y = sinc(x-(D));
hold on 
plot(x,sinc_y,LineWidth=lw);
xlabel('Samples',Interpreter="latex",FontSize=fontSize),
ylabel('Impulse response',Interpreter="latex",FontSize=fontSize)
xlim([-1 6])
grid

title('N = 5',Interpreter="latex",FontSize=fontSize)


subplot(3,2,3)
stem(n_opt, h_int_opt,LineWidth=lw),
x = linspace(0,6,100);
sinc_y = sinc(x-(D));
hold on 
plot(x,sinc_y,LineWidth=lw);
xlabel('Samples',Interpreter="latex",FontSize=fontSize),
ylabel('Impulse response',Interpreter="latex",FontSize=fontSize)
xlim([-1 7])
grid
title('N = 6',Interpreter="latex",FontSize=fontSize)


subplot(3,2,4)
stem(n_7, h_int_7,LineWidth=lw),
x = linspace(0,7,100);
sinc_y = sinc(x-(D));
hold on 
plot(x,sinc_y,LineWidth=lw);
xlabel('Samples',Interpreter="latex",FontSize=fontSize),
ylabel('Impulse response',Interpreter="latex",FontSize=fontSize)
xlim([-1 8])
grid
title('N = 7',Interpreter="latex",FontSize=fontSize)


subplot(3,2,5)
stem(n_8, h_int_8,LineWidth=lw),
x = linspace(0,8,100);
sinc_y = sinc(x-(D));
hold on 
plot(x,sinc_y,LineWidth=lw);
xlabel('Samples',Interpreter="latex",FontSize=fontSize),
ylabel('Impulse response',Interpreter="latex",FontSize=fontSize)
title('Lagrange interpolation',Interpreter="latex",FontSize=fontSize)
xlim([-1 9])
grid
title('N = 8',Interpreter="latex",FontSize=fontSize)


subplot(3,2,6)
stem(n_9, h_int_9,LineWidth=lw),
x = linspace(0,9,100);
sinc_y = sinc(x-(D));
hold on
plot(x,sinc_y,LineWidth=lw);
xlabel('Samples',Interpreter="latex",FontSize=fontSize),
ylabel('Impulse response',Interpreter="latex",FontSize=fontSize)
title('N = 9',Interpreter="latex",FontSize=fontSize)
xlim([-1 10])
grid
sgtitle('Impulse response Comparison',Interpreter="latex",FontSize=fontSize+3);
saveas(gcf,"Plots/FiltersIR.png")
%[H_int, f] = freqz(h_int_1, 1, 'half', fs);

figure('Renderer', 'painters', 'Position', [10 10 1000 1200])

subplot(3,1,1)
subtitle('Frequency response',Interpreter="latex",FontSize=fontSize);
% plot(f_1*fs/pi/2, abs(H_int_1));
% hold on
% plot(f_2*fs/pi/2, abs(H_int_2));
% hold on
%norm_freq = (0,1,f_4*fs/pi/2);
plot(f_4/(pi), abs(H_int_4),LineWidth=lw);
hold on
plot(f_5/(pi), abs(H_int_5),LineWidth=lw);
hold on
plot(f_opt/(pi), abs(H_int_opt),LineWidth=lw);
hold on
plot(f_7/(pi), abs(H_int_7),LineWidth=lw);
hold on
plot(f_8/(pi), abs(H_int_8),LineWidth=lw);
hold on
plot(f_9/(pi), abs(H_int_9),LineWidth=lw);
ylim([0 1.7]);
%xlim([0 fs/2]);
xlabel('Frequency [Hz]',Interpreter="latex",FontSize=fontSize), ylabel('Magnitude [dB]',Interpreter="latex",FontSize=fontSize);
%legend(["N=1","N=2","N=3","N= floor(2*D)","N=7","N= floor(4*D)","N= floor(5*D)"],'Location','northwest',Interpreter="latex");
legend(["N=4","N=5","N=6","N=7","N=8","N=9"],'Location','northwest',Interpreter="latex",FontSize=fontSize);
hold on
%axis([0, fs/2, 0 2]),
grid

subplot(3,1,2)
subtitle('Group Delay',Interpreter="latex",FontSize=fontSize);
% grpdelay(h_int_1, 1, 1024, 'half', fs);
% hold on 
% grpdelay(h_int_2, 1, 1024, 'half', fs);
% hold on 
[gd,w]=grpdelay(h_int_4, 1, 1024);
plot(w/pi,gd,LineWidth=lw)

hold on
[gd,w]=grpdelay(h_int_5, 1, 1024);
plot(w/pi,gd,LineWidth=lw)

hold on 
[gd,w]=grpdelay(h_int_opt, 1, 1024);
plot(w/pi,gd,LineWidth=lw)

hold on 
[gd,w]=grpdelay(h_int_7, 1, 1024);
plot(w/pi,gd,LineWidth=lw)

hold on 
[gd,w]=grpdelay(h_int_8, 1, 1024);
plot(w/pi,gd,LineWidth=lw)

hold on 
[gd,w]=grpdelay(h_int_9, 1, 1024);
plot(w/pi,gd,LineWidth=lw)

ylim([1 4]);
xlabel('Normalized Frequency $(x \pi [rad/sample])$',Interpreter="latex",FontSize=fontSize), ylabel('Group Delay [samples]',Interpreter="latex",FontSize=fontSize);
%xlim([0 fs/2]);
%legend(["N=1","N=2","N=3","N= floor(2*D)","N=7","N= floor(4*D)","N= floor(5*D)"],'Location','northwest',Interpreter="latex");
legend(["N=4","N=5","N=6","N=7","N=8","N=9"],'Location','northwest',Interpreter="latex",FontSize=fontSize);
grid
subplot(3,1,3)
subtitle('Phase Delay',Interpreter="latex",FontSize=fontSize);
[pd,w]=phasedelay(h_int_4, 1, 1024);
plot(w/pi,pd,LineWidth=lw)

hold on 
[pd,w]=phasedelay(h_int_5, 1, 1024);
plot(w/pi,pd,LineWidth=lw)

hold on 
[pd,w]=phasedelay(h_int_opt, 1, 1024);
plot(w/pi,pd,LineWidth=lw)

hold on 
[pd,w]=phasedelay(h_int_7, 1, 1024);
plot(w/pi,pd,LineWidth=lw)

hold on 
[pd,w]=phasedelay(h_int_8, 1, 1024);
plot(w/pi,pd,LineWidth=lw)

hold on 
[pd,w]=phasedelay(h_int_9, 1, 1024);
plot(w/pi,pd,LineWidth=lw)
xlabel('Normalized Frequency $(x \pi [rad/sample])$',Interpreter="latex",FontSize=fontSize), ylabel('Phase Delay [samples]',Interpreter="latex",FontSize=fontSize);
grid
%xlim([0 1024]);
%legend(["N=1","N=2","N=3","N= floor(2*D)","N=7","N= floor(4*D)","N= floor(5*D)"],'Location','northwest',Interpreter="latex");
legend(["N=4","N=5","N=6","N=7","N=8","N=9"],'Location','northwest',Interpreter="latex",FontSize=fontSize);
sgtitle('Filter response orders comparison',Interpreter="latex",FontSize=fontSize);
saveas(gcf,"Plots/FreqResponseComparison.png")


x = linspace(0,2*floor(D),100);
sinc_y = sinc(x-(D));

figure('Renderer', 'painters', 'Position', [10 10 1000 800])
stem(n_opt, h_int_opt),
hold on 
plot(x,sinc_y);