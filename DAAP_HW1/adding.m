% =========================================================================
% adding function for restoring audio signal from short-time
% windowType gives the overlap assigned in the windowing function
% 
%
% =========================================================================



function [addedSignal] = adding(shortTimeSignal,windowType,windowLength)

M = windowLength;
if(strcmp(windowType,"rectwin"))
     R = M;
elseif(strcmp(windowType,"bartlett"))
    R = M/2;
elseif(strcmp(windowType,"hamming"))
    R = M/2;
elseif(strcmp(windowType,"hann"))
    R= M/2;
else 
    disp("wrong name use of rect window!")
    R=M;
end
[windowLength,chunks_num] = size(shortTimeSignal);
signalLength = (chunks_num * R)+(M-R);
addedSignal = zeros(signalLength,1);
ii=1;
for nn = 1:R:(signalLength-(M-R))
    s_chunk = shortTimeSignal(:,ii);
    addedSignal( nn:nn+(M-1),1) = addedSignal(nn:nn+(M-1),1) + s_chunk;
    ii =ii+1;
    %disp(strcat("nn=", num2str(nn)))
end

if(strcmp(windowType,"rectwin"))
elseif(strcmp(windowType,"bartlett"))
    signal_rem = zeros(length(addedSignal)-M,1);
    signal_rem = signal_rem + addedSignal((M/2+1):(M/2+length(signal_rem)),1);
    addedSignal = signal_rem;
elseif(strcmp(windowType,"hamming"))
    signal_rem = zeros(length(addedSignal)-M,1);
    signal_rem = signal_rem + addedSignal((M/2+1):(M/2+length(signal_rem)),1);
    addedSignal = signal_rem;
elseif(strcmp(windowType,"hann"))
    signal_rem = zeros(length(addedSignal)-M,1);
    signal_rem = signal_rem + addedSignal((M/2+1):(M/2+length(signal_rem)),1);
    addedSignal = signal_rem;
else
end

%audiowrite("addedSignal.wav",addedSignal,44100);