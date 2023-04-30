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