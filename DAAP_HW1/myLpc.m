function [shapingFilters,whiteningFilters] =  myLpc(st_signal,taps_number,steepDesc,numberOfCycles)

if(steepDesc==1)
    p = taps_number;
    [M,N] = size(st_signal); % M window length, N number of windows
    r = zeros(p+1, N);
    
    for nn = 1:N
        for ii_kk = 0:p
            for mm = 1:M-(ii_kk)
                r(ii_kk+1,nn) = r(ii_kk+1,nn)+ st_signal(mm,nn)*st_signal(mm+ii_kk,nn);
            end
        end 
    end
    
    r_0 = r(1:end-1,:);
    r_1 = r(2:end,:);
    R = zeros(p,p,N);
    a = zeros(p, N);
    grad_J = zeros(p,N);
    for nn = 1:N
        R(:,:,nn) = toeplitz(r_0(:,nn));
         % adding new method for steepest descent, I think it's wrong
         [~,lambda] = eig(R(:,:,nn));
         mu = 2/(max(diag(lambda))*10);
        % for ss = 1:numberOfCycles
            grad_J(:,nn) = -(2 *r_1(:,nn)) + (2* R(:,:,nn) * a(:,nn)); 
            a(:,nn) = a(:,nn) - (0.5 * mu *grad_J(:,nn));
%         while norm(grad_J) > sqrt(p*0.1)
%             grad_J(:,nn) = -2 *r_1(:,nn) + (2* R(:,:,nn) * a(:,nn)); 
% 
%             a(:,nn+1) = a(:,nn) + (0.5 * mu * - grad_J(1,nn));
%         end
    end
    
    for nn = 1:N
        while norm(grad_J) > sqrt(p*100)
            R_a = (2* R(:,:,nn) * a(:,nn)); 
            grad_J(:,nn) = -(2 *r_1(:,nn)) + R_a;

            a(:,nn) = a(:,nn) - (0.5 * mu * ( grad_J(:,nn)));
            
            %a(:,nn) = a(:,nn) + ( mu * ( grad_J(:,nn)));
            
            norm(grad_J)
        end
    end
    a_1 = zeros(p+1,N);
    A = zeros(M,N);
    H = zeros(M,N);
    
    for nn = 1: N
        a_1(:,nn) = vertcat(1, -a(:,nn));
        
        [A(:,nn),~] = freqz(a_1(:,nn),1,"whole",M);
        %A(:,nn) = A(:,nn)/(mean(A(:,nn)/mean(st_signal(:,nn))));
        [H(:,nn),w] = freqz(1,a_1(:,nn),"whole",M);
        %H(:,nn) = H(:,nn)/(mean(H(:,nn))/mean(st_signal(:,nn)));
        %A(:,nn) = 1/(H(:,nn));
        
    end 
    shapingFilters = H;
    
    whiteningFilters = A;

else
    p = taps_number;
    [M,N] = size(st_signal); % M window length, N number of windows
    r = zeros(p+1, N);
    
    for nn = 1:N
        for ii_kk = 0:p
            for mm = 1:M-(ii_kk)
                r(ii_kk+1,nn) = r(ii_kk+1,nn)+ st_signal(mm,nn)*st_signal(mm+ii_kk,nn);
            end
        end 
    end
    
    r_0 = r(1:end-1,:);
    r_1 = r(2:end,:);
    R = zeros(p,p,N);
    a = zeros(p, N);
    for nn = 1:N
        R(:,:,nn) = toeplitz(r_0(:,nn));
        a(:,nn) = R(:,:,nn)\r_1(:,nn);
    end
    
    a_1 = zeros(p+1,N);
    A = zeros(M,N);
    H = zeros(M,N);
    
    for nn = 1: N
        a_1(:,nn) = vertcat(1, -a(:,nn));
        
    [A(:,nn),~] = freqz(a_1(:,nn),1,"whole",M);
        %A(:,nn) = A(:,nn)/(mean(A(:,nn)/mean(st_signal(:,nn))));
        [H(:,nn),w] = freqz(1,a_1(:,nn),"whole",M);
        %H(:,nn) = H(:,nn)/(mean(H(:,nn))/mean(st_signal(:,nn)));
        %A(:,nn) = 1/(H(:,nn));
        
    end 
    shapingFilters = H;
    
    whiteningFilters = A;
end