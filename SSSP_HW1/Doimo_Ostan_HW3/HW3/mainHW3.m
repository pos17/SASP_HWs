%------------------------------------------%
%        *** SSSP - HOMEWORK #3 ***        %
%------------------------------------------%
%     Emulation of the Leslie Speaker      %
%------------------------------------------%
% Name:       Alberto Doimo                %
% Student ID: 10865196                     %
% Name:       Paolo Ostan                  %
% Student ID: 10868276                     %
%------------------------------------------%

clear; close all; clc;

%% modulation speed
mod_speed = 'tremolo';

%% Read the input file
[x, Fs] = audioread('HammondRef.wav');
x=x(:,1);           % take left channel only

%% FX parameters 
switch lower(mod_speed)
    
    case {'chorale'}
        freq=2;
            
    case {'tremolo'}
        freq=6;
        
    otherwise
        error('mod_speed \"%s\" not found.', mod_speed)

end

%% Apply FX
[y,y_lpf,y_hpf,y_hp_sdf]= leslie(x, Fs, freq);

%figure 
%plot(y)
%% Avoid any (possible) clipping
y = rescale(y,-1.,1.);

%% Playback
%audiowrite([mod_speed,'.wav'], y, Fs);
soundsc(y, Fs)

%% Read the reference audio file
dir_name = 'Leslie_ref';
addpath(dir_name);
[y_ref, ~] = audioread(fullfile(dir_name, strcat(mod_speed,'.wav')));

%% Display the MSE
MSE = mean(abs(y-y_ref).^2);
MSE_str = sprintf('MSE: %g', MSE);
disp(MSE_str)
