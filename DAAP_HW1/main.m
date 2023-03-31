function [] = main(instrFileName,speechFileName, outputFileName, ...
    solMode,tuningMu,minThresh,cycNumMax,initialValues,resample_var,verbose,plotCom)

%=========================================================================%
%                           DAAP HW1 main                                 %
%
%   solMode = "steepDesc" or "linSolve"                                                                         %  
%=========================================================================%


%% xcorr, circshift, cell, window

if not(isfolder("audioOutputs"))
    mkdir("audioOutputs")
end

addpath('audioInputs')
addpath('audioOutputs')

[instr_t,instr_Fs] = audioread(instrFileName);
[speech_t,speech_Fs] = audioread(speechFileName);

%% Resample
if(resample_var == 1)
    fsin    = instr_Fs;
    fsout   = 16000;
    m       = lcm(fsin,fsout);
    up      = m/fsin;
    down    = m/fsout;
    instr_t    = resample(instr_t, up, down);
    instr_Fs=fsout;
    
    fsin    = speech_Fs;
    fsout   = 16000;
    m       = lcm(fsin,fsout);
    up      = m/fsin;
    down    = m/fsout;
    speech_t    = resample(speech_t, up, down);
    speech_Fs=fsout;
end
%%

speech_t(end+1:length(instr_t),1) =0;
instr_t_len =length(instr_t);

taps_speech = 44;
taps_music = 80;

wl =  2048; 

[instr_st_signal,chunksNum_instr] = windowing(instr_t,"hann",wl,instr_Fs,verbose);
[speech_st_signal,chunksNum_speech] = windowing(speech_t,"hann",wl,speech_Fs,verbose);

if strcmp(solMode,"steepDesc") 
    [instr_H,instr_A] =  myLpc(instr_st_signal,taps_music,"steepDesc",tuningMu,minThresh,cycNumMax,0,initialValues,verbose);
    [speech_H,speech_A] =  myLpc(speech_st_signal,taps_speech,"steepDesc",tuningMu,minThresh,cycNumMax,1,initialValues,verbose);
elseif strcmp(solMode,"linSolve") 
    [instr_H,instr_A] =  myLpc(instr_st_signal,taps_music,"linSolve",tuningMu,minThresh,cycNumMax,0,initialValues,verbose);
    [speech_H,speech_A] =  myLpc(speech_st_signal,taps_speech,"linSolve",tuningMu,minThresh,cycNumMax,0,initialValues,verbose);
elseif strcmp(solMode,"mixed1") 
    [instr_H,instr_A] =  myLpc(instr_st_signal,taps_music,"linSolve",tuningMu,minThresh,cycNumMax,0,initialValues,verbose);
    [speech_H,speech_A] =  myLpc(speech_st_signal,taps_speech,"steepDesc",tuningMu,minThresh,cycNumMax,1,initialValues,verbose);
end

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
% M_instr_H = max(abs(instr_H), [], 'all');
% M_instr_sig_st = max(abs(instr_st_signal_w), [], 'all');
% M_speech_H = max(abs(speech_H), [], 'all');
% M_speech_sig_st = max(abs(speech_st_signal_w), [], 'all');
% 
% instr_H = (instr_H/M_instr_H)*M_instr_sig_st;
% speech_H = (speech_H/M_speech_H)*M_speech_sig_st;    



instr_st_res =  zeros(wl, chunksNum_instr); 
instr_st_res_w = zeros(wl,chunksNum_instr);
talking_instr_st_res =  zeros(wl, chunksNum_instr); 
talking_instr_st_res_w = zeros(wl,chunksNum_instr);
for nn = 1:chunksNum_instr
    %instr_A(:,nn) = instr_A(:,nn)/(mean(instr_A(:,nn)/mean(instr_st_signal_w(:,nn))));
    
    %instr_st_res_w(:,nn) =  instr_st_signal_w(:,nn).* instr_A(:,nn);
    instr_st_res_w(:,nn) =  instr_st_signal_w(:,nn)./ instr_H(:,nn);
    
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

talking_instr_lin = adding(talking_instr_st_res,"hann",wl);
instr_lin = adding(instr_st_res,"hann",wl);
%st_res_lin = adding(instr_st_res,0.5,wl);
%st_res_lin = adding(instr_st_res,0.5,wl);
%st_res_lin = st_res_lin / mean(st_res_lin);
%st_res_lin = adding(instr_st_res,0.5,wl);


%% testing speech shaping filter 

% sine = linspace(1,length(speech_t)/speech_Fs,length(speech_t)).';
% sine_shaped = zeros(length(speech_t),1);
% sine = 0.1 *sin(500*sine);
% sine_shaped_st = zeros(wl,chunksNum_speech);
% sine_shaped_st_w = zeros(wl,chunksNum_speech);
% sine_st_w = zeros(wl,chunksNum_speech);
% sine_st = windowing(sine,"hamming",wl);
% for nn = 1:chunksNum_speech
%     sine_st_w(:,nn) = fft(sine_st(:,nn));
%     sine_shaped_st_w(:,nn) = sine_st_w(:,nn) .* speech_H(:,nn);
%     sine_shaped_st(:,nn) = ifft(sine_shaped_st_w(:,nn));
% end
% 
% sine_shaped = adding(sine_shaped_st,0.5,wl);
% sine_shaped = sine_shaped/max(abs(sine_shaped));
% audiowrite("./audioOutputs/"+"whiteShaped.wav",real(sine_shaped),instr_Fs);

%%

talking_instr_lin =talking_instr_lin /max(abs(talking_instr_lin));
audiowrite("./audioOutputs/"+outputFileName,real(talking_instr_lin),instr_Fs);
%audiowrite("./audioOutputs/"+instr_res_name,real(instr_lin),instr_Fs);


%% PLOTS

% variables initialization

axlabelsize = 15;
titlesize = 22;
legendsize = 15;

% Instrument vs filter spectrum comparison 

w = linspace(0,instr_Fs/2,wl/2);
figure('Renderer', 'painters', 'Position', [10 10 1000 600])

for i=1:9
    subplot(3,3,i); 
    index=i*floor(chunksNum_instr/9);
    plot(w,10*log10(abs(instr_st_signal_w(1:end/2,index)).^2));
    xlabel("frequency [Hz]",Interpreter="latex");
    ylabel("[dB]",Interpreter="latex");
    hold on 
    plot(w,10*log10(abs(instr_H(1:end/2,index)).^2),"Linestyle","-", "Linewidth",1.3)
    title(index,Interpreter='Latex');
    legend("signal chunk","filter shape",Interpreter='Latex');
    %legend('boxoff');
end
sgtitle('Instrument filter comparison, ' + plotCom, FontSize=titlesize, Interpreter='Latex');

% Speech vs filter spectrum comparison 

w = linspace(0,speech_Fs/2,wl/2);
figure('Renderer', 'painters', 'Position', [10 10 1000 600])

for i=1:9
    subplot(3,3,i); 
    index=i*floor(chunksNum_instr/9);
    plot(w,10*log10(abs(speech_st_signal_w(1:end/2, index)).^2));
    xlabel("frequency [Hz]",Interpreter="latex");
    ylabel("[dB]",Interpreter="latex");
    hold on 
    plot(w,10*log10(abs(speech_H(1:end/2,index).^2)),"Linestyle","-", "Linewidth",1.3);
    title(index,Interpreter='Latex');
    legend("signal chunk","filter shape",Interpreter='Latex');
    %legend('boxoff');

end
sgtitle('Speech filter comparison'+plotCom, FontSize=titlesize, Interpreter='Latex');

% Power spectrograms

figure('Renderer', 'painters', 'Position', [10 10 1000 700])

subplot(311);
w1 = linspace(0,instr_Fs/2,wl/2);
t =linspace(0,instr_t_len/instr_Fs,chunksNum_instr);
instr_st_signal_w_ps = instr_st_signal_w(1:wl/2,:);
surf(t,w1,(abs(instr_st_signal_w_ps).^2)./(max(abs(instr_st_signal_w_ps).^2)),EdgeColor="none");
ylim([0 2000]);
view(0,90);
xlabel("time [s]",Interpreter="latex");
ylabel("frequency [Hz]",Interpreter="latex");
title('Piano',Interpreter='Latex');

subplot(312)
w1 = linspace(0,speech_Fs/2,wl/2);
t =linspace(0,instr_t_len/instr_Fs,chunksNum_speech);
speech_st_signal_w_ps = speech_st_signal_w(1:wl/2,:);
surf(t,w1,(abs(speech_st_signal_w_ps).^2)./(max(abs(speech_st_signal_w_ps).^2)),EdgeColor="none");
%ylim([0 3000])
view(0,90);
xlabel("time [s]",Interpreter="latex");
ylabel("frequency [Hz]",Interpreter="latex");
title('Speech',Interpreter='Latex');
 
subplot(313)
w1 = linspace(0,speech_Fs/2,wl/2);
t =linspace(0,instr_t_len/instr_Fs,chunksNum_speech);
talking_instr_st_res_w_ps = talking_instr_st_res_w(1:wl/2,:);
surf(t,w1,abs(talking_instr_st_res_w_ps).^2./(max(abs(talking_instr_st_res_w_ps).^2)),EdgeColor="none");
%ylim([0 3000])
view(0,90);
xlabel("time [s]",Interpreter="latex");
ylabel("frequency [Hz]",Interpreter="latex");
title('Talking instrument',Interpreter='Latex');

sgtitle('Spectrogram', FontSize=titlesize, Interpreter='Latex');

end