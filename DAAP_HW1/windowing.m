function [shortTimeSignal,numOfChunks] = windowing(signal,windowType,windowLength,Fs,verbose)

%==========================================================================

%                   function for windowing the signal
%
%==========================================================================

nZeroPad =0;
M = windowLength;
if(strcmp(windowType,"rectwin"))
    w = rectwin(windowLength);
    R = M;
elseif(strcmp(windowType,"bartlett"))
    w = bartlett(windowLength);
    R = M/2;
    % zero padding adding half windowlenght    
    % before signal samples to have full signal
    signal_pad = zeros(M + length(signal),1);
    signal_pad(M/2+1:M/2+length(signal)) = signal;

elseif(strcmp(windowType,"hamming"))
    w = hamming(windowLength,"periodic");
    R = M/2;
    signal_pad = zeros(M + length(signal),1);
    signal_pad(M/2+1:M/2+length(signal)) = signal;
elseif(strcmp(windowType,"hann"))
    w = hann(windowLength,"periodic");
    R= M/2;
    signal_pad = zeros(M + length(signal),1);
    %sg = signal_pad(M/2+1:M/2+length(signal));
    signal_pad(M/2+1:M/2+length(signal)) = signal;
    signal = signal_pad;
else 
    disp("wrong name use of rect window!")
    w = rectwin(windowLength);
    R=M;
end
signalLength = length(signal);
% cola conditions need to be respected!

%t_buckets = ceil(signalLength/windowLength);

signal(end+1:end+1+M) = 0;
chunks_num = ceil(signalLength / R);

shortTimeSignal = zeros(M,chunks_num);

ii=1;
for nn = 1:R:signalLength
    s_chunk = signal(nn:nn+(M-1));
    s_chunk1 = s_chunk.* w;
    shortTimeSignal(:,ii) = s_chunk1;
    ii =ii+1;
end
if (verbose == 1) 
    shortTimeWindows = zeros(M,chunks_num);    
    for nn = 1:chunks_num
        shortTimeWindows(:,nn) = w;
    end
    % sum test before adding 
    addedWindows = adding(shortTimeWindows,windowType,windowLength);
    figure
    plot(linspace(0,length(addedWindows)/Fs,length(addedWindows)),addedWindows,"b-o");
    title("TEST FOR COLA CONDITION")
    ylim([0.9,1.1]);
    xlim([-0.3,length(addedWindows)/Fs+0.3]);
end

% figure
% plot(w)
% hold on
% plot(shortTimeSignal(:,200))
numOfChunks = chunks_num;