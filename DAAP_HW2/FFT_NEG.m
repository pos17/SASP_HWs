function X = FFT_NEG(Y) 
    n = length(Y);
    if(n==1)
        X = Y;
    else 
        y_0 = Y(1:2:length(Y));
        y_1 = Y(2:2:length(Y));
        a_0 = myifft(y_0);
        a_1 = myifft(y_1);
        ON_NEG = exp(-1i*2*pi/n);
        ON_NEGI = 1;
        a = zeros(1,n);
        for ii =1:((n/2))
            t = ON_NEGI * a_1(ii);
            a(1,ii) = a_0(ii) + t;
            a(ii+n/2) = a_0(ii) - t;
            ON_NEGI = ON_NEGI * ON_NEG; 
        end 
        X = a;
    end

end