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
y_1_t_len = length(y_1_t);
y_2_t_len = length(y_2_t);
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



%test stft 
M = 1024;
R = M/2;
Y_1_st = mystft(y_1_t,"hann",M,R,M);
winsNum = size(Y_1_st,2);

figure 
w1 = linspace(0,y_1_Fs/2,M/2);
t =linspace(0,y_1_t_len/y_1_Fs,winsNum);
Y_11_st = Y_1_st(1:M/2,:);
surf(t,w1,abs(Y_11_st)./(max(abs(Y_11_st))),EdgeColor="none");
%ylim([0 2000])
view(0,90);
xlabel("time [s]",Interpreter="latex");
ylabel("frequency [Hz]",Interpreter="latex");
title('spectrogram',Interpreter='Latex');

%test istft

y_1_st_1_res = myistft(Y_1_st,2*M,R);

audiowrite("./audioOutputs/"+"mystftTest.wav",real(y_1_st_1_res)/max(real(y_1_st_1_res)),y_1_Fs);

figure
plot(real(y_1_st_1_res))

% algorithm 

y_1_st = mystft(y_1_)



