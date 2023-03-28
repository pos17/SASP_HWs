clc; close all; clear all;

%==========================================================================
%                           DAAP HW1 main
%                           
%==========================================================================

% xcorr, circshift, cell, window
addpath('audioInputs')
addpath('audioOutputs')
instr_name = "piano.wav";
instr_name_recon = "piano_recon.wav";
speech_name = "speech.wav";
talking_instr_res_name = "talking_instr_res.wav";
instr_res_name = "instr_filter_res.wav";
st_res_sub_name = "subtracted_res.wav";
st_res_sub_name_test = "subtracted_res.wav";


[instr_t,instr_Fs] = audioread(instr_name);
[speech_t,speech_Fs] = audioread(speech_name);

speech_t(end+1:length(instr_t),1) =0;
instr_t_len =length(instr_t);

taps_speech = 800;
taps_music = 800;

wl =  1024; 

[instr_st_signal,chunksNum_instr] = windowing(instr_t,"hann",wl);
[speech_st_signal,chunksNum_speech] = windowing(speech_t,"hann",wl);


[instr_H,instr_A] =  myLpc(instr_st_signal,taps_music);
[speech_H,speech_A] =  myLpc(speech_st_signal,taps_speech);


instr_st_signal_w = zeros(wl,chunksNum_instr);
speech_st_signal_w = zeros(wl,chunksNum_speech);

for nn = 1:chunksNum_instr
    instr_st_signal_w(:,nn) = fft(instr_st_signal(:,nn));
    speech_st_signal_w(:,nn) = fft(speech_st_signal(:,nn));
    
end


instr_st_res =  zeros(wl, chunksNum_instr); 
instr_st_res_w = zeros(wl,chunksNum_instr);
talking_instr_st_res =  zeros(wl, chunksNum_instr); 
talking_instr_st_res_w = zeros(wl,chunksNum_instr);
for nn = 1:chunksNum_instr
    %instr_A(:,nn) = instr_A(:,nn)/(mean(instr_A(:,nn)/mean(instr_st_signal_w(:,nn))));
    instr_st_res_w(:,nn) =  instr_st_signal_w(:,nn).* instr_A(:,nn);
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

audiowrite("./audioOutputs/"+talking_instr_res_name,real(talking_instr_lin),instr_Fs);
audiowrite("./audioOutputs/"+instr_res_name,real(instr_lin),instr_Fs);
