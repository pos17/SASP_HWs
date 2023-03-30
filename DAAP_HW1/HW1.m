%==========================================================================
%                           DAAP HW1
%                  Talking Instrument creation                         
%==========================================================================

clc; close all; clear all;

instr_name = "piano.wav";
speech_name = "speech.wav";
talking_instr_res_name = "talking_instr_res.wav";
instr_res_name = "instr_filter_res.wav";
st_res_sub_name = "subtracted_res.wav";
st_res_sub_name_test = "subtracted_res.wav";

avviso = audioread("audioInputs/avviso.mp3");
%f = @() main(instr_name, speech_name, talking_instr_res_name,"linSolve"); % handle to function
%timeit(f,1)
disp("time for speech steepest descent from previous filter")
tic
main(instr_name, speech_name, talking_instr_res_name,"mixed1",0.3,0.5*10^-1,100000,"prev",0);
%main(instr_name, speech_name, talking_instr_res_name,"linSolve");
toc
disp("time for speech steepest descent equals to zero")

tic
main(instr_name, speech_name, talking_instr_res_name,"mixed1",0.3,0.5*10^-1,100000,"zeros",0);
%main(instr_name, speech_name, talking_instr_res_name,"linSolve");
toc
disp("time for speech steepest descent equals to one")

tic
main(instr_name, speech_name, talking_instr_res_name,"mixed1",0.3,0.5*10^-1,100000,"ones",0);
%main(instr_name, speech_name, talking_instr_res_name,"linSolve");
toc

disp("time for speech steepest descent equals to one")


tic
main(instr_name, speech_name, talking_instr_res_name,"mixed1",0.3,0.5*10^-1,100000,"randn",0);
%main(instr_name, speech_name, talking_instr_res_name,"linSolve");
toc

sound(avviso(1:5*44100),44100);