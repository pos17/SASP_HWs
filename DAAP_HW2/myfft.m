function Y = myfft(X) 
    n = length(X);
    if(n==1)
        Y = X;
    else
        a_0 = X(1:2:length(X));
        a_1 = X(2:2:length(X));
        y_0 = myfft(a_0);
        y_1 = myfft(a_1);

        ON = exp((2*pi*1i)/n);
        ONI = 1;
        y = zeros(1,n);
        for ii =1:((n/2))
            t = ONI * y_1(ii);
            y(1,ii) = y_0(ii) + t;
            y(ii+n/2) = y_0(ii) - t;
            ONI = ONI * ON; 
        end 
        Y = y;
    end
end