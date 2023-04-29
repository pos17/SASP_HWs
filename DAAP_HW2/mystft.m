function [Y_st] = mystft(x,windowType,windowLength,overlap,zeropadding)
    P=zeropadding;
    M =windowLength;
    x = x(:);
    R=overlap;
    Nframes = floor((length(x)-M)/R) + 1;
    %X = zeros(N/2+1,Nframes);
    Y_st = zeros(M+P,Nframes);
    x(length(x)+1:(Nframes*R)+(2*M-R))=0;
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

    for n = 0 : Nframes
        xn = x(n*R+1 : n*R+M);         % Extract a block of the signal
        yn = w .* xn;
        yn(end+1:end+P) = 0; % Window the signal block
        temp = fft(yn);              % Compute the DFT of the windowed block
        Y_st(:, n+1) = temp;  %temp(1 : N/2+1) ;   % Store the value of the DFT in the output matrix  
    end


end