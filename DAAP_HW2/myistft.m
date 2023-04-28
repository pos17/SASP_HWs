function [y] = myistft(x,windowLength,overlap)
    R=overlap;
    M=windowLength;
    if(windowLength ~= size(x,1) )
        disp("something wrong!!!!!!!");
    end
    Nframes = size(x,2);
    y_t_length = (Nframes*R)+(windowLength-R);
    y_t = zeros(y_t_length,1);
    

    for m = 0 : Nframes-1
        temp = ifft(x(:,m+1));  
        y_t(m*R+1 : m*R+M)=y_t(m*R+1 : m*R+M) + temp;
    end
    y=y_t;

end