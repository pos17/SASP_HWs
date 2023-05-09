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
N = 1;

h_int_1 = LagrangeInterpolation(N,D);
n_1 = 0:N;


[H_int_1, f_1] = freqz(h_int_1, 1, 'half', fs);


N = 2;

h_int_2 = LagrangeInterpolation(N,D);
n_2 = 0:N;
[H_int_2, f_2] = freqz(h_int_2, 1, 'half', fs);

N = 3;

h_int_3 = LagrangeInterpolation(N,D);
n_3= 0:N;
[H_int_3, f_3] = freqz(h_int_3, 1, 'half', fs);

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
N = 4*floor(D);
n_2opt = 0:N;
h_int_2opt = LagrangeInterpolation(N,D);
[H_int_2opt, f_2opt] = freqz(h_int_2opt, 1, 'half', fs);

% optimal filter settings
N = 5*floor(D);
n_3opt = 0:N;
h_int_3opt = LagrangeInterpolation(N,D);
[H_int_3opt, f_3opt] = freqz(h_int_3opt, 1, 'half', fs);



%% Lagrange interpolation plots


figure('Renderer', 'painters', 'Position', [10 10 1000 800])
stem(n_1, h_int_1),
hold on 
stem(n_2, h_int_2),
hold on 
stem(n_3, h_int_3),
hold on 
stem(n_opt, h_int_opt),
hold on 
stem(n_7, h_int_7),
hold on 
stem(n_2opt, h_int_2opt),
hold on 
stem(n_3opt, h_int_3opt),
 
xlabel('Samples',Interpreter="latex"),
ylabel('Impulse response',Interpreter="latex")
title('Lagrange interpolation',Interpreter="latex")
legend(["N=1","N=2","N=3","N= floor(2*D)","N=7","N= floor(4*D)","N= floor(5*D)"],'Location','northwest',Interpreter="latex");
saveas(gcf,"Plots/FiltersIR.png")
%[H_int, f] = freqz(h_int_1, 1, 'half', fs);

figure('Renderer', 'painters', 'Position', [10 10 1000 800])
subplot(2,1,1)
plot(f_1*fs/pi/2, abs(H_int_1));
hold on
plot(f_2*fs/pi/2, abs(H_int_2));
hold on
plot(f_3*fs/pi/2, abs(H_int_3));
hold on
plot(f_opt*fs/pi/2, abs(H_int_opt));
hold on
plot(f_7*fs/pi/2, abs(H_int_7));
hold on
plot(f_2opt*fs/pi/2, abs(H_int_2opt));
hold on
plot(f_3opt*fs/pi/2, abs(H_int_3opt));
xlabel('Frequency [Hz]',Interpreter="latex"), ylabel('Magnitude [dB]',Interpreter="latex");
legend("N=1","N=2","N=3","N= floor(2*D)","N=7","N= floor(4*D)","N= floor(5*D)",'Location','northwest',Interpreter="latex");
hold on
axis([0, fs/2, 0 2]), grid

subplot(2,1,2)
grpdelay(h_int_1, 1, 1024, 'half', fs);
hold on 
grpdelay(h_int_2, 1, 1024, 'half', fs);
hold on 
grpdelay(h_int_3, 1, 1024, 'half', fs);
hold on 
grpdelay(h_int_opt, 1, 1024, 'half', fs);
hold on 
grpdelay(h_int_7, 1, 1024, 'half', fs);
hold on 
grpdelay(h_int_2opt, 1, 1024, 'half', fs);
hold on 
grpdelay(h_int_3opt, 1, 1024, 'half', fs);
legend("N=1","N=2","N=3","N= floor(2*D)","N=7","N= floor(4*D)","N= floor(5*D)",'Location','northwest',Interpreter="latex");
sgtitle('Frequency response for orders comparison',Interpreter="latex");
saveas(gcf,"Plots/FreqResponseComparison.png")