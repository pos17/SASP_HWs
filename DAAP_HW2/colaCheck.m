signal = ones(80000,1);

M = 512;
R = M/4;
P = M+1;
winType = "hann";
sig_stft = mystft(signal,winType,M,R,P);

sig_istft = myistft(sig_stft,2*M+1,R);

figure 

plot(sig_istft)
ylim([0.98 1.02])



% fft check

if not(isfolder("audioOutputs"))
    mkdir("audioOutputs")
end

addpath('audioInputs')
addpath('audioOutputs')

y_1_name = "y1.wav";
[y_1_t,y_1_Fs] = audioread(y_1_name);

y_1_t = y_1_t(1:4096);


y_1_w = fft(y_1_t);

y_1_w_ang= angle(y_1_w);

y_1_w_ang_1 = y_1_w_ang(1:2048,1);
y_1_w_ang_2 = flip(y_1_w_ang(2049:end),1);

figure
subplot(2,1,1)
plot(abs(y_1_w))
subplot(2,1,2)
plot(angle(y_1_w))

figure
plot(y_1_w_ang_1)
hold on
plot(-y_1_w_ang_2)
