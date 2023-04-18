function [Y_st] = mystft(x,windowType,windowLength,overlap)
    N =windowLength;
    R = overlap;
    x = x(:);
    R=overlap;
    Nframes = floor((length(x)-N)/R) + 1;
    %X = zeros(N/2+1,Nframes);
    Y_st = zeros(N,Nframes);
    x(length(x)+1:(Nframes*R)+(2*N-R))=0;
    %windowing
    if(strcmp(windowType,"rectwin"))
        w = rectwin(windowLength);
    elseif(strcmp(windowType,"bartlett"))
        w = bartlett(windowLength);
    elseif(strcmp(windowType,"hamming"))
        w = hamming(windowLength,"periodic");
    elseif(strcmp(windowType,"hann"))
        w = hann(windowLength,"periodic");
    else 
        disp("wrong name use of rect window!")
        w = rectwin(windowLength);
    end

    w = w(:);

    for m = 0 : Nframes
        xm = x(m*R+1 : m*R+N);         % Extract a block of the signal
        ym = w .* xm;                   % Window the signal block
        temp = fft(ym);              % Compute the DFT of the windowed block
        Y_st(:, m+1) = temp;  %temp(1 : N/2+1) ;   % Store the value of the DFT in the output matrix  
    end


end