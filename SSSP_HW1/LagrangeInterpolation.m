function [h_int] = LagrangeInterpolation(N,D)
    h_int = zeros(N+1,1);
    for nn = 0:N
        k = (0:N)';
        k = k(k~=nn);
        h_int(nn+1) = prod((D-k)./(nn-k));
    end
end 