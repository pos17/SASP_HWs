%% DFT, CONVOLUTION and STFT - EXERCISE 8

function [X, Nframes] = STFT(x, w, R, N)
% This function performs the Short-Time Fourier Transform of signal x
%   x - signal (usually very long)
%   w - frame window
%   R - frame overlap
%   N - frame length
% It returns the STFT Matrix X and the number of frames Nframes

    % Make both the signal and the window column vectors
    x = x(:);
    w = w(:);

    % Compute number of frames in x
    Nframes = floor((length(x)-N)/R) + 1 ;
    
    % Window length
    M = length(w);

    % Initialize STFT Matrix
    X = zeros(N/2+1,Nframes);   % Matrix size: (#frequency bins, #frames)
                                % Notice that #frequency bins = N/2 + 1
                                % -> we store the one-sided spectrum

    % Compute the STFT
    for m = 0 : Nframes
        xm = x(m*R+1 : m*R+M) ;         % Extract a block of the signal
        ym = w .* xm;                   % Window the signal block
        temp = fft(ym, N);              % Compute the DFT of the windowed block
        X(:, m+1) = temp(1 : N/2+1) ;   % Store the value of the DFT in the output matrix
    end

end

