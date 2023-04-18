clc
close all 
clear all
% Try STFT on a signal loaded from a .wav file

% Load the signal
[y_1, Fs_1] = audioread('audiofiles/y1.wav');  % built-in function audioread
                                   % audioread returns the signal x and the
                                   % corresponding sampling frequency Fs
y_1 = y_1(:, 1);                       % Takes only the left channel (in order to have a mono signal) 

% STFT Parameters
M = floor(0.05 * Fs_1);   % Window length (50msec)
                         % Notice that M is [#samples] -> M = duration * Fs
R = floor(M/4);          % Hop size (overlap 25%)
                         % M and R must be integers -> use floor for
                         % rounding in a "conservative" way
N = 2^(ceil(log2(M)));   % (nearest power of 2) > M -> we can use also nextpow2
w = hanning(M)';         % Hanning window of length M

% STFT computation
[y_1, N_1] = STFT(y_1, w, R, N);

% Plot the SpectroGram
t = (0 : N_1 - 1) * R / Fs_1 ;    % Time axis [sec]
                                    % Take into account that frames are overlapped 
f = Fs_1 * (0 : N_1/2) / N_1;             % Frequency axis [Hz]
                                    % Take into account that the spectrum
                                    % is one-sided -> from 0 to Fs/2, with step Fs/N
figure
imagesc(t, f, 10*log10(abs(y_1).^2)); % dB scale
                                    % Notice that we are plotting |X|^2 (power of the signal)
colorbar; 
%caxis([-150 50]);
xlabel('time (sec)');
ylabel('frequency (Hz)');
axis xy
%In MATLAB exists a built-in function called spectrogram that computes the STFT
% Built-in Spectrogram
% [Xs, fax, tax] = spectrogram(y_1, w, M-R+1, N_1, Fs_1);   % Inputs:
                                                    % - Signal
                                                    % - Window
                                                    % - #frames
                                                    % - #FFT bins
                                                    % - Sampling Frequency                                          
% imagesc(t, f, 10*log10(abs(Xs).^2));
% xlabel('time (sec)');
% ylabel('frequency (Hz)');
% axis xy

%%
clear
[y_2, Fs_2] = audioread('audiofiles/y2.wav');
y_2 = y_2(:, 1);  
% STFT Parameter
M = floor(0.05 * Fs_2);    
R = floor(M/4);       
N = 2^(ceil(log2(M)));   
w = hanning(M)';     
[y_2, N_2] = STFT(y_2, w, R, N);

% Plot the SpectroGram
t = (0 : N_2 - 1) * R / Fs_2 ;    % Time axis [sec]
                                    % Take into account that frames are overlapped 
f = Fs_2 * (0 : N_2/2) / N_2;             % Frequency axis [Hz]
                                    % Take into account that the spectrum
                                    % is one-sided -> from 0 to Fs/2, with step Fs/N
figure
imagesc(t, f, 10*log10(abs(y_2).^2)); % dB scale
                                    % Notice that we are plotting |X|^2 (power of the signal)
colorbar;
%caxis([-150 50]);
xlabel('time (sec)');
ylabel('frequency (Hz)');
axis xy
%In MATLAB exists a built-in function called spectrogram that computes the STFT
% Built-in Spectrogram
%[Xs, fax, tax] = spectrogram(y_2, w, M-R+1, N, Fs_2);   % Inputs:
                                                    % - Signal
                                                    % - Window
                                                    % - #frames
                                                    % - #FFT bins
                                                    % - Sampling Frequency                                          
% imagesc(t, f, 10*log10(abs(Xs).^2));
% xlabel('time (sec)');
% ylabel('frequency (Hz)');
% axis xy