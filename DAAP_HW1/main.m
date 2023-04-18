function [] = main(instrFileName,speechFileName, outputFileName,...
    wl_speech,wl_instr,taps_speech,taps_music,solMode,tuningMu,...
    minThresh,cycNumMax,initialValues,resample_var,spec,verbose,plotCom,...
    normalize,plotName)
    warning('off');
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

%taps_speech = 44;
%taps_music = 80;

%wl_speech =  2048; 

[instr_st_signal,chunksNum_instr] = windowing(instr_t,"hann",wl_instr,instr_Fs,verbose,strcat("instr_",plotName));
[speech_st_signal,chunksNum_speech] = windowing(speech_t,"hann",wl_speech,speech_Fs,verbose,strcat("speech_",plotName));

if strcmp(solMode,"steepDesc") 
    if(verbose == 1)
       disp("Precision check plots for piano signal")
    end
    [instr_H,instr_A] =  myLpc(instr_st_signal,taps_music,"steepDesc",tuningMu,minThresh,cycNumMax,0,initialValues,instr_Fs,verbose,plotName);
    if(verbose == 1)
       disp("Precision check plots for speech signal")
    end
    [speech_H,speech_A] =  myLpc(speech_st_signal,taps_speech,"steepDesc",tuningMu,minThresh,cycNumMax,1,initialValues,speech_Fs,verbose,plotName);
elseif strcmp(solMode,"linSolve") 
    if(verbose == 1)
       disp("Precision check plots for piano signal")
    end
    [instr_H,instr_A] =  myLpc(instr_st_signal,taps_music,"linSolve",tuningMu,minThresh,cycNumMax,0,initialValues,instr_Fs,verbose,plotName);
    if(verbose == 1)
       disp("Precision check plots for speech signal")
    end
    [speech_H,speech_A] =  myLpc(speech_st_signal,taps_speech,"linSolve",tuningMu,minThresh,cycNumMax,0,initialValues,speech_Fs,verbose,plotName);
elseif strcmp(solMode,"mixed1") 
    if(verbose == 1)
       disp("Precision check plots for piano signal")
    end
    [instr_H,instr_A] =  myLpc(instr_st_signal,taps_music,"linSolve",tuningMu,minThresh,cycNumMax,0,initialValues,instr_Fs,verbose,plotName);
    if(verbose == 1)
       disp("Precision check plots for speech signal")
    end
    [speech_H,speech_A] =  myLpc(speech_st_signal,taps_speech,"steepDesc",tuningMu,minThresh,cycNumMax,1,initialValues,speech_Fs,verbose,plotName);
end

instr_st_signal_w = zeros(wl_instr,chunksNum_instr);
instr_st_res = zeros(wl_instr,chunksNum_instr);
instr_st_res_w = zeros(wl_instr,chunksNum_instr);
speech_st_signal_w = zeros(wl_speech,chunksNum_speech);
normValuesInstr = ones(1,chunksNum_instr);
for nn = 1:chunksNum_instr
    instr_st_signal_w(:,nn) = fft(instr_st_signal(:,nn));

    if(normalize == "mean")
        normValue = mean(abs(instr_H(:,nn)))/mean(abs(instr_st_signal_w(:,nn)));
        %instr_H(:,nn) = (instr_H(:,nn)/mean(abs(instr_H(:,nn))))*mean(abs(instr_st_signal_w(:,nn)));
        %normValuesInstr(1,nn) = mean(abs(instr_H(:,nn)))/mean(abs(instr_st_signal_w(:,nn)));
        instr_H(:,nn) = instr_H(:,nn)/normValue;
        normValuesInstr(1,nn) = normValue;
    elseif(normalize == "max")
        %instr_H(:,nn) = (instr_H(:,nn)/max(abs(instr_H(:,nn))))*max(abs(instr_st_signal_w(:,nn)));
        %normValuesInstr(1,nn) = max(abs(instr_H(:,nn)))/max(abs(instr_st_signal_w(:,nn)));
        normValue = max(abs(instr_H(:,nn)))/max(abs(instr_st_signal_w(:,nn)));
        instr_H(:,nn) = instr_H(:,nn)/normValue;
        normValuesInstr(1,nn) = normValue;
    
    elseif(normalize == "none")
    end

    instr_st_res_w(:,nn) =  instr_st_signal_w(:,nn)./ instr_H(:,nn);   
    instr_st_res(:,nn) = ifft(instr_st_res_w(:,nn));
   
    
end

instr_lin = adding(instr_st_res,"hann",wl_instr);
instr_lin = instr_lin(1:length(speech_t));
[instr_st_signal_speech,chunksNum_instr_speech] = windowing(instr_lin,"hann",wl_speech,instr_Fs,verbose,strcat("instr_",plotName));

%instr_st_res =  zeros(wl_instr, chunksNum_instr); 
%instr_st_res_w = zeros(wl_instr,chunksNum_instr);
talking_instr_st_res =  zeros(wl_speech, chunksNum_instr); 
talking_instr_st_res_w = zeros(wl_speech,chunksNum_instr);
normValuesSpeech = ones(1,chunksNum_instr_speech);
for nn = 1:chunksNum_instr_speech
   speech_st_signal_w(:,nn) = fft(speech_st_signal(:,nn));
   instr_st_signal_speech_w(:,nn) = fft(instr_st_signal_speech(:,nn));
   if(normalize == "mean")
        normValue = mean(abs(speech_H(:,nn)))/mean(abs(speech_st_signal_w(:,nn)));
        %speech_H(:,nn) = (speech_H(:,nn)/mean(abs(speech_H(:,nn))))*mean(abs(speech_st_signal_w(:,nn)));
        speech_H(:,nn) = speech_H(:,nn)/normValue;
        normValuesSpeech(1,nn) = normValue;
   elseif(normalize == "max")
        %speech_H(:,nn) = (speech_H(:,nn)/max(abs(speech_H(:,nn))))*max(abs(speech_st_signal_w(:,nn)));   
        normValue = max(abs(speech_H(:,nn)))/max(abs(speech_st_signal_w(:,nn)));
        speech_H(:,nn) = speech_H(:,nn)/normValue;
        normValuesSpeech(1,nn) = normValue;
   
        %normValuesSpeech(1,nn) = max(abs(speech_H(:,nn)))/max(abs(speech_st_signal_w(:,nn)));
   elseif(normalize == "none")
   end
   talking_instr_st_res_w(:,nn) = instr_st_signal_speech_w(:,nn) .* speech_H(:,nn);
   talking_instr_st_res(:,nn) = ifft(talking_instr_st_res_w(:,nn));  
end




% M_instr_H = max(abs(instr_H), [], 'all');
% M_instr_sig_st = max(abs(instr_st_signal_w), [], 'all');
% M_speech_H = max(abs(speech_H), [], 'all');
% M_speech_sig_st = max(abs(speech_st_signal_w), [], 'all');
% 
% instr_H = (instr_H/M_instr_H)*M_instr_sig_st;
% speech_H = (speech_H/M_speech_H)*M_speech_sig_st;    




% for nn = 1:chunksNum_instr
%     instr_st_res_w(:,nn) =  instr_st_signal_w(:,nn)./ instr_H(:,nn);   
%     talking_instr_st_res_w(:,nn) = instr_st_res_w(:,nn) .* speech_H(:,nn);
%     instr_st_res(:,nn) = ifft(instr_st_res_w(:,nn));
%     talking_instr_st_res(:,nn) = ifft(talking_instr_st_res_w(:,nn));
% end

% try to normalize output

%for nn = 1:chunksNum_instr
%    instr_st_res(:,nn) = instr_st_res(:,nn) / mean(instr_st_res(:,nn));
%end


%st_res_lin = reshape(st_res,[t_buckets *wl 1]);

talking_instr_lin = adding(talking_instr_s t_res,"hann",wl_speech);

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

% Steepest descend comparison

w = linspace(0,speech_Fs/2,wl_speech/2);
figure('Renderer', 'painters', 'Position', [10 10 1000 600])

index=ceil(chunksNum_speech/2);
plot(w,10*log10(abs(speech_st_signal_w(1:end/2, index)).^2));
xlabel("frequency [Hz]",Interpreter="latex");
ylabel("[dB]",Interpreter="latex");
xlim([0 speech_Fs/2]);
hold on 
plot(w,10*log10(abs(speech_H(1:end/2,index).^2)),"Linestyle","-", "Linewidth",1.3);
legend("signal chunk","filter shape",Interpreter='Latex');
grid minor
%legend('boxoff');

sgtitle('Speech filter comparison' +plotCom +index, FontSize=titlesize, Interpreter='Latex');
plotName = 'Speech filter comparison' +plotCom +index;
saveas(gcf,strcat("plots/",plotName),"png");
% Instrument vs filter spectrum comparison 

w = linspace(0,instr_Fs/2,wl_instr/2);
figure('Renderer', 'painters', 'Position', [10 10 1000 600])
disp("Instr");
for i=1:9
    subplot(3,3,i); 
    index=i*floor(chunksNum_instr/9);
    plot(w,10*log10(abs(instr_st_signal_w(1:end/2,index)).^2));
    xlabel("frequency [Hz]",Interpreter="latex");
    ylabel("[dB]",Interpreter="latex");
    xlim([0 instr_Fs/2]);
    hold on 
    plot(w,10*log10(abs(instr_H(1:end/2,index)).^2),"Linestyle","-", "Linewidth",1.3)
    title(index,Interpreter='Latex');
    legend("signal chunk","filter shape",Interpreter='Latex');
    grid minor
    %legend('boxoff');
    disp(strcat("plot: ", num2str(index)));
    disp(strcat("normalization value: ", num2str(normValuesInstr(1,index))));
end
sgtitle('Instrument filter comparison' + plotCom, FontSize=titlesize, Interpreter='Latex');
saveas(gcf,strcat("plots/InstrFilterComp_",plotName),"png");
% Speech vs filter spectrum comparison 

w = linspace(0,speech_Fs/2,wl_speech/2);
figure('Renderer', 'painters', 'Position', [10 10 1000 600])
disp("Speech");
for i=1:9
    subplot(3,3,i); 
    index=i*floor(chunksNum_speech/9);
    plot(w,10*log10(abs(speech_st_signal_w(1:end/2, index)).^2));
    xlabel("frequency [Hz]",Interpreter="latex");
    ylabel("[dB]",Interpreter="latex");
    xlim([0 speech_Fs/2]);
    hold on 
    plot(w,10*log10(abs(speech_H(1:end/2,index).^2)),"Linestyle","-", "Linewidth",1.3);
    title(index,Interpreter='Latex');
    legend("signal chunk","filter shape",Interpreter='Latex');
    grid minor
    %legend('boxoff');
    disp(strcat("plot: ", num2str(index)));
    disp(strcat("normalization value: ", num2str(normValuesSpeech(1,index))));
end
sgtitle('Speech filter comparison' +plotCom, FontSize=titlesize, Interpreter='Latex');
saveas(gcf,strcat("plots/SpeechFilterComp_",plotName),"png");

% Power spectrograms
if spec == 1
figure('Renderer', 'painters', 'Position', [10 10 1000 700])

subplot(311);
w1 = linspace(0,instr_Fs/2,wl_instr/2);
t =linspace(0,instr_t_len/instr_Fs,chunksNum_instr);
instr_st_signal_w_ps = instr_st_signal_w(1:wl/2,:);
surf(t,w1,(abs(instr_st_signal_w_ps).^2)./(max(abs(instr_st_signal_w_ps).^2)),EdgeColor="none");
ylim([0 2000]);
view(0,90);
xlabel("time [s]",Interpreter="latex");
ylabel("frequency [Hz]",Interpreter="latex");
title('Piano',Interpreter='Latex');

subplot(312)
w1 = linspace(0,speech_Fs/2,wl_speech/2);
t =linspace(0,instr_t_len/instr_Fs,chunksNum_speech);
speech_st_signal_w_ps = speech_st_signal_w(1:wl/2,:);
surf(t,w1,(abs(speech_st_signal_w_ps).^2)./(max(abs(speech_st_signal_w_ps).^2)),EdgeColor="none");
ylim([0 2000])
view(0,90);
xlabel("time [s]",Interpreter="latex");
ylabel("frequency [Hz]",Interpreter="latex");
title('Speech',Interpreter='Latex');
 
subplot(313)
w1 = linspace(0,speech_Fs/2,wl_speech/2);
t =linspace(0,instr_t_len/instr_Fs,chunksNum_speech);
talking_instr_st_res_w_ps = talking_instr_st_res_w(1:wl/2,:);
surf(t,w1,abs(talking_instr_st_res_w_ps).^2./(max(abs(talking_instr_st_res_w_ps).^2)),EdgeColor="none");
ylim([0 2000])
view(0,90);
xlabel("time [s]",Interpreter="latex");
ylabel("frequency [Hz]",Interpreter="latex");
title('Talking instrument',Interpreter='Latex');

sgtitle('Spectrogram', FontSize=titlesize, Interpreter='Latex');
end
end