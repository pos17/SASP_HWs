clc; close all; clear all;
axlabelsize = 15;
titlesize = 22;
legendsize = 15;
%==========================================================================
%                           DAAP HW1 main
%                           
%==========================================================================

% xcorr, circshift, cell, window

instr_name = "piano.wav";
instr_name_recon = "piano_recon.wav";
speech_name = "speech.wav";
talking_instr_res_name = "talking_instr_res.wav";
instr_res_name = "instr_filter_res.wav";
st_res_sub_name = "subtracted_res.wav";
st_res_sub_name_test = "subtracted_res.wav";

if not(isfolder("audioOutputs"))
    mkdir("audioOutputs")
end

addpath('audioInputs')
addpath('audioOutputs')

[instr_t,instr_Fs] = audioread(instr_name);
[speech_t,speech_Fs] = audioread(speech_name);

%% Resample

fsin    = instr_Fs;
fsout   = 16000;
m       = lcm(fsin,fsout);
up      = m/fsin;
down    = m/fsout;
instr_t    = resample(instr_t, up, down);
%audiowrite([a_filename,'_22050','.wav'], x_22, fsout);
instr_Fs=fsout;

fsin    = speech_Fs;
fsout   = 16000;
m       = lcm(fsin,fsout);
up      = m/fsin;
down    = m/fsout;
speech_t    = resample(speech_t, up, down);
speech_Fs=fsout;

%%

speech_t(end+1:length(instr_t),1) =0;
instr_t_len =length(instr_t);

taps_speech = 100;
taps_music = 30;

wl =  1024; 

[instr_st_signal,chunksNum_instr] = windowing(instr_t,"hann",wl);
[speech_st_signal,chunksNum_speech] = windowing(speech_t,"hann",wl);


[instr_H,instr_A] =  myLpc(instr_st_signal,taps_music,0,10);
[speech_H,speech_A] =  myLpc(speech_st_signal,taps_speech,0,10);


instr_st_signal_w = zeros(wl,chunksNum_instr);
speech_st_signal_w = zeros(wl,chunksNum_speech);

for nn = 1:chunksNum_instr
    instr_st_signal_w(:,nn) = fft(instr_st_signal(:,nn));
    speech_st_signal_w(:,nn) = fft(speech_st_signal(:,nn));
    
end


for nn = 1:chunksNum_instr
    instr_H(:,nn) = (instr_H(:,nn)/max(abs(instr_H(:,nn))))*max(abs(instr_st_signal_w(:,nn)));
    speech_H(:,nn) = (speech_H(:,nn)/max(abs(speech_H(:,nn))))*max(abs(speech_st_signal_w(:,nn)));
    
end





instr_st_res =  zeros(wl, chunksNum_instr); 
instr_st_res_w = zeros(wl,chunksNum_instr);
talking_instr_st_res =  zeros(wl, chunksNum_instr); 
talking_instr_st_res_w = zeros(wl,chunksNum_instr);
for nn = 1:chunksNum_instr
    %instr_A(:,nn) = instr_A(:,nn)/(mean(instr_A(:,nn)/mean(instr_st_signal_w(:,nn))));
    
    instr_st_res_w(:,nn) =  instr_st_signal_w(:,nn).* instr_A(:,nn);
    %instr_st_res_w(:,nn) =  instr_st_signal_w(:,nn)./ instr_H(:,nn);
    
    talking_instr_st_res_w(:,nn) = instr_st_res_w(:,nn) .* speech_H(:,nn);
    %st_res_w(:,nn) = st_signal_w(:,nn)./ H(:,nn);
    instr_st_res(:,nn) = ifft(instr_st_res_w(:,nn));
    talking_instr_st_res(:,nn) = ifft(talking_instr_st_res_w(:,nn));
end

% try to normalize output

%for nn = 1:chunksNum_instr
%    instr_st_res(:,nn) = instr_st_res(:,nn) / mean(instr_st_res(:,nn));
%end


%st_res_lin = reshape(st_res,[t_buckets *wl 1]);

talking_instr_lin = adding(talking_instr_st_res,0.5,wl);
instr_lin = adding(instr_st_res,0.5,wl);
%st_res_lin = adding(instr_st_res,0.5,wl);
%st_res_lin = adding(instr_st_res,0.5,wl);
%st_res_lin = st_res_lin / mean(st_res_lin);
%st_res_lin = adding(instr_st_res,0.5,wl);


%% testing speech shaping filter 

sine = linspace(1,length(speech_t)/speech_Fs,length(speech_t)).';
sine_shaped = zeros(length(speech_t),1);
sine = 0.1 *sin(500*sine);
sine_shaped_st = zeros(wl,chunksNum_speech);
sine_shaped_st_w = zeros(wl,chunksNum_speech);
sine_st_w = zeros(wl,chunksNum_speech);
sine_st = windowing(sine,"hamming",wl);
for nn = 1:chunksNum_speech
    sine_st_w(:,nn) = fft(sine_st(:,nn));
    sine_shaped_st_w(:,nn) = sine_st_w(:,nn) .* speech_H(:,nn);
    sine_shaped_st(:,nn) = ifft(sine_shaped_st_w(:,nn));
end

sine_shaped = adding(sine_shaped_st,0.5,wl);
sine_shaped = sine_shaped/max(abs(sine_shaped));
audiowrite("./audioOutputs/"+"whiteShaped.wav",real(sine_shaped),instr_Fs);

%%

talking_instr_lin =talking_instr_lin /max(abs(talking_instr_lin));
audiowrite("./audioOutputs/"+talking_instr_res_name,real(talking_instr_lin),instr_Fs);
audiowrite("./audioOutputs/"+instr_res_name,real(instr_lin),instr_Fs);


%% plots

%for nn = 1:chunksNum_instr
%    instr_H(:,nn) = (instr_H(:,nn) /max(instr_H(:,nn)))*max(instr_st_signal_w(:,nn));
%end

w = linspace(instr_Fs/2,-instr_Fs/2,wl);
figure('Renderer', 'painters', 'Position', [10 10 1000 600])

subplot(3,3,1); 
plot(w,10*log10(abs(instr_st_signal_w(:,50).^2)))
%xlim([0,instr_Fs/2])
xlabel(" []")
ylabel(" []")
hold on 
plot(w,10*log10(abs(instr_H(:,50)).^2),"LineStyle","--")
title("50")
legend("signal chunk","filter shape","test filter shape")


subplot(3,3,2); 
plot(w,10*log10(abs(instr_st_signal_w(:,100)).^2))
xlabel(" []")
ylabel(" []")
hold on 
plot(w,10*log10(abs(instr_H(:,100)).^2),"LineStyle","--")
title("100")
legend("signal chunk","filter shape","test filter shape")

subplot(3,3,3); 
plot(w,10*log10(abs(instr_st_signal_w(:,150)).^2))
xlabel(" []")
ylabel(" []")
hold on 
plot(w,10*log10(abs(instr_H(:,150)).^2),"LineStyle","--")
title("150")
legend("signal chunk","filter shape","test filter shape")

subplot(3,3,4); 
plot(w,10*log10(abs(instr_st_signal_w(:,200)).^2))
xlabel(" []")
ylabel(" []")
hold on 
plot(w,10*log10(abs(instr_H(:,200)).^2),"LineStyle","--")
title("200")
legend("signal chunk","filter shape","test filter shape")

subplot(3,3,5); 
plot(w,10*log10(abs(instr_st_signal_w(:,250)).^2))
xlabel(" []")
ylabel(" []")
hold on 
plot(w,10*log10(abs(instr_H(:,250)).^2),"LineStyle","--")
title("250")
legend("signal chunk","filter shape","test filter shape")

subplot(3,3,6); 
plot(w,10*log10(abs(instr_st_signal_w(:,300)).^2))
xlabel(" []")
ylabel(" []")
hold on 
plot(w,10*log10(abs(instr_H(:,300)).^2),"LineStyle","--")
title("300")
legend("signal chunk","filter shape","test filter shape")

subplot(3,3,7); 
plot(w,10*log10(abs(instr_st_signal_w(:,350)).^2))
xlabel(" []")
ylabel(" []")
hold on 
plot(w,10*log10(abs(instr_H(:,350)).^2),"LineStyle","--")
title("350")
legend("signal chunk","filter shape","test filter shape")

subplot(3,3,8); 
plot(w,10*log10(abs(instr_st_signal_w(:,400)).^2))
xlabel(" []")
ylabel(" []")
hold on 
plot(w,10*log10(abs(instr_H(:,400)).^2),"LineStyle","--")
title("400")
legend("signal chunk","filter shape","test filter shape")

subplot(3,3,9); 
plot(w,10*log10(abs(instr_st_signal_w(:,450)).^2))
xlabel(" []")
ylabel(" []")
hold on 
plot(w,10*log10(abs(instr_H(:,450)).^2),"LineStyle","--")
title("450")
legend("signal chunk","filter shape","test filter shape")

sgtitle('Instrument filter comparison', FontSize=titlesize)



figure
title("speech filter comparison")
plot(w,10*log10(abs(speech_st_signal_w(:,100)).^2))
%plot(w,abs(st_signal_w(:,600)))
hold on 
plot(w,10*log10(abs(speech_H(:,100)).^2),"LineStyle","--")
%plot(w,abs(H(:,600)))
hold on 
%plot(w,10*log10(abs(H_test(:,600)).^2),"LineStyle","--")
%plot(w,abs(H(:,600)))



w = linspace(-speech_Fs/2,speech_Fs/2,wl);
figure('Renderer', 'painters', 'Position', [10 10 1000 600])

subplot(3,3,1); 
plot(w,10*log10(abs(speech_st_signal_w(:,50)).^2))
xlabel(" []")
ylabel(" []")
hold on 
plot(w,10*log10(abs(speech_H(:,50)).^2),"LineStyle","--")
title("50")
legend("signal chunk","filter shape","test filter shape")

subplot(3,3,2); 
plot(w,10*log10(abs(speech_st_signal_w(:,100)).^2))
xlabel(" []")
ylabel(" []")
hold on 
plot(w,10*log10(abs(speech_H(:,100)).^2),"LineStyle","--")
title("100")
legend("signal chunk","filter shape","test filter shape")

subplot(3,3,3); 
plot(w,10*log10(abs(speech_st_signal_w(:,150)).^2))
xlabel(" []")
ylabel(" []")
hold on 
plot(w,10*log10(abs(speech_H(:,150)).^2),"LineStyle","--")
title("150")
legend("signal chunk","filter shape","test filter shape")

subplot(3,3,4); 
plot(w,10*log10(abs(speech_st_signal_w(:,200)).^2))
xlabel(" []")
ylabel(" []")
hold on 
plot(w,10*log10(abs(speech_H(:,200)).^2),"LineStyle","--")
title("200")
legend("signal chunk","filter shape","test filter shape")

subplot(3,3,5); 
plot(w,10*log10(abs(speech_st_signal_w(:,250)).^2))
xlabel(" []")
ylabel(" []")
hold on 
plot(w,10*log10(abs(speech_H(:,250)).^2),"LineStyle","--")
title("250")
legend("signal chunk","filter shape","test filter shape")

subplot(3,3,6); 
plot(w,10*log10(abs(speech_st_signal_w(:,300)).^2))
xlabel(" []")
ylabel(" []")
hold on 
plot(w,10*log10(abs(speech_H(:,300)).^2),"LineStyle","--")
title("300")
legend("signal chunk","filter shape","test filter shape")

subplot(3,3,7); 
plot(w,10*log10(abs(speech_st_signal_w(:,350)).^2))
xlabel(" []")
ylabel(" []")
hold on 
plot(w,10*log10(abs(speech_H(:,350)).^2),"LineStyle","--")
title("350")
legend("signal chunk","filter shape","test filter shape")

subplot(3,3,8); 
plot(w,10*log10(abs(speech_st_signal_w(:,400)).^2))
xlabel(" []")
ylabel(" []")
hold on 
plot(w,10*log10(abs(speech_H(:,400)).^2),"LineStyle","--")
title("400")
legend("signal chunk","filter shape","test filter shape")

subplot(3,3,9); 
plot(w,10*log10(abs(speech_st_signal_w(:,450)).^2))
xlabel(" []")
ylabel(" []")
hold on 
plot(w,10*log10(abs(speech_H(:,450)).^2),"LineStyle","--")
title("450")
legend("signal chunk","filter shape","test filter shape")

sgtitle('speech filter comparison', FontSize=titlesize)




%xlim([0,1])
legend("signal chunk","filter shape","test filter shape")


figure 
plot(talking_instr_lin)