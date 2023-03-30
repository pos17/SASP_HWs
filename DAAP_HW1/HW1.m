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

main(instr_name, speech_name, talking_instr_res_name,"steepDesc");