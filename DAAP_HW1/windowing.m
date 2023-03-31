function [shortTimeSignal,numOfChunks] = windowing(signal,windowType,windowLength,verbose)

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
elseif(strcmp(windowType,"hamming"))
    w = hamming(windowLength,"periodic");
    R = M/2;
elseif(strcmp(windowType,"hann"))
    w = hann(windowLength,"periodic");
    R= M/2;
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
        shortTimeWindows(:,ii) = w;
    end
    % sum test before adding 
    addedWindows = adding(shortTimeWindows,R,windowLength);
    figure
    plot(1:length(addedWindows),addedWindows);
    title("TEST TO SEE IF COLA CONDITION IS RESPECTED")
end

% figure
% plot(w)
% hold on
% plot(shortTimeSignal(:,200))
numOfChunks = chunks_num;