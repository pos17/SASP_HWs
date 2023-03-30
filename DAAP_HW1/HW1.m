%==========================================================================
%                           DAAP HW1
%                  Talking Instrument creation                         
%==========================================================================

clc; close all; clear all;
axlabelsize = 15;
titlesize = 22;
legendsize = 15;




instr_name = "piano.wav";
speech_name = "speech.wav";
talking_instr_res_name = "talking_instr_res.wav";
instr_res_name = "instr_filter_res.wav";
st_res_sub_name = "subtracted_res.wav";
st_res_sub_name_test = "subtracted_res.wav";

avviso = audioread("audioInputs/avviso.mp3");
%f = @() main(instr_name, speech_name, talking_instr_res_name,"linSolve"); % handle to function
%timeit(f,1)
%disp("time for speech steepest descent from previous filter")
disp("threshold 10^-1")
tic
main(instr_name, speech_name, talking_instr_res_name,"mixed1",0.3,0.5*10^-1,100000,"zeros",1,0);
%main(instr_name, speech_name, talking_instr_res_name,"linSolve");
toc
%disp("time for speech steepest descent equals to zero")

%disp("time for speech steepest descent from previous filter")
disp("threshold 10^-2")
tic
main(instr_name, speech_name, talking_instr_res_name,"mixed1",0.3,0.5*10^-2,100000,"zeros",1,0);
%main(instr_name, speech_name, talking_instr_res_name,"linSolve");
toc
disp("threshold 10^-4")

%disp("time for speech steepest descent from previous filter")
tic
main(instr_name, speech_name, talking_instr_res_name,"mixed1",0.3,0.5*10^-4,100000,"zeros",1,0);
%main(instr_name, speech_name, talking_instr_res_name,"linSolve");
toc
disp("time for speech steepest descent equals to zero")




sound(avviso(1:5*44100),44100);
%%
tic
main(instr_name, speech_name, talking_instr_res_name,"mixed1",0.3,0.5*10^-1,100000,"zeros",1,0);
%main(instr_name, speech_name, talking_instr_res_name,"linSolve");
toc
disp("time for speech steepest descent equals to one")

tic
main(instr_name, speech_name, talking_instr_res_name,"mixed1",0.3,0.5*10^-1,100000,"ones",0);
%main(instr_name, speech_name, talking_instr_res_name,"linSolve");
toc

disp("time for speech steepest descent equals to random values")


tic
main(instr_name, speech_name, talking_instr_res_name,"mixed1",0.3,0.5*10^-1,100000,"randn",0);
%main(instr_name, speech_name, talking_instr_res_name,"linSolve");
toc

sound(avviso(1:5*44100),44100);