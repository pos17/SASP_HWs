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


%% FIRST TEST WITH VERBOSE PARAMETERS TO CHECK FUNCTIONING OF THE SCRIPT
tic
main(instr_name, speech_name, talking_instr_res_name,"linSolve",0.3,10^-1,100000,"zeros",1,1,"threshold 10^(-1)");
toc

%% test for various threshold values using steepest descent
close all;clc;
disp("threshold 10^-1")
tic
main(instr_name, speech_name, talking_instr_res_name,"mixed1",0.3,10^-1,100000,"zeros",1,0,"threshold 10^(-1)");
toc
disp("threshold 10^-2")
tic
main(instr_name, speech_name, talking_instr_res_name,"mixed1",0.3,10^-2,100000,"zeros",1,0,"threshold 10^(-2)");
toc
disp("threshold 10^-4")

tic
main(instr_name, speech_name, talking_instr_res_name,"mixed1",0.3,10^-4,100000,"zeros",1,0,"threshold 10^(-4)");
toc



sound(avviso(1:5*44100),44100);

%% test for various coefficients initialization values using steepest descent
tic
main(instr_name, speech_name, talking_instr_res_name,"mixed1",0.3,0.5*10^-1,100000,"zeros",1,0);
toc

tic
main(instr_name, speech_name, talking_instr_res_name,"mixed1",0.3,0.5*10^-1,100000,"zeros",1,0);
toc
disp("time for speech steepest descent equals to one")

tic
main(instr_name, speech_name, talking_instr_res_name,"mixed1",0.3,0.5*10^-1,100000,"ones",0);
toc

disp("time for speech steepest descent equals to random values")


tic
main(instr_name, speech_name, talking_instr_res_name,"mixed1",0.3,0.5*10^-1,100000,"randn",0);
toc

sound(avviso(1:5*44100),44100);