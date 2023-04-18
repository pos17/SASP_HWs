clc; close all; clear all;

y_1_name = "y1.wav";
y_2_name = "y2.wav";

if not(isfolder("audioOutputs"))
    mkdir("audioOutputs")
end

addpath('audioInputs')
addpath('audioOutputs')

[y_1_t,y_1_Fs] = audioread(y_1_name);
[y_2_t,y_2_Fs] = audioread(y_2_name);

pw = ceil(log(length(y_1_t))/log(2));
y_1_t(length(y_1_t):2^pw)=0;

Y_1_w = myfft(y_1_t);
figure 

plot(10*log10(abs(Y_1_w).^2));
hold on 
plot(10*log10(abs(fft(y_1_t)).^2),"--");

y_1_w_ref = fft(y_1_t);
y_1_w_ref = y_1_w_ref';
y_1_sub = abs(Y_1_w) - abs(y_1_w_ref);
figure

plot(y_1_sub);

y_1_t_res = myifft(Y_1_w);
figure
plot(real(y_1_t_res));

audiowrite("./audioOutputs/"+"myfftTest.wav",real(y_1_t_res)/max(real(y_1_t_res)),y_1_Fs);

y_1_t_res_real = real(y_1_t_res)/max(real(y_1_t_res));
y_1_t_real = real(y_1_t)/max(real(y_1_t));

y_1_t_sub =y_1_t_res_real - y_1_t_real';

figure

plot(y_1_t_sub);

y_1_t_sub = y_1_t_sub / max(y_1_t_sub);
audiowrite("./audioOutputs/"+"myfftTestSUB.wav",real(y_1_t_sub)/max(real(y_1_t_sub)),y_1_Fs);


