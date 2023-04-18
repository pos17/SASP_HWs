function x = myifft(Y)
    n = length(Y);
    a = FFT_NEG(Y);
    a = a/n;
    x=a;
end