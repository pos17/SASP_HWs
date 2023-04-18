%% xcorr, circshift, cell, window

if not(isfolder("audioOutputs"))
    mkdir("audioOutputs")
end

addpath('audioInputs')
addpath('audioOutputs')

[instr_t,instr_Fs] = audioread(instrFileName);
[speech_t,speech_Fs] = audioread(speechFileName);