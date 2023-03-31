function [shapingFilters,whiteningFilters] =  myLpc(st_signal,taps_number,...
 solveMode,tuningMu,minThresh,cycNumMax,convergenceTest,initialValues,verbose)
axlabelsize = 15;
titlesize = 22;
legendsize = 15;



% initialValues:    ones array of ones
%                   zeros array of zeros
%                   prev previous values of the previous chunk
if strcmp(solveMode,"steepDesc") 
   
   comparing_conv_1 = 1;
   comparing_conv_2 = 2;
    p = taps_number;
    [M,N] = size(st_signal); % M window length, N number of windows
    chosen_conv = ceil(N/2);
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
    if(strcmp(initialValues,"zeros"))
        a = zeros(p, N);
    elseif (strcmp(initialValues,"ones"))
        a = ones(p, N);
    elseif(strcmp(initialValues,"prev"))
        a = zeros(p, N);
    elseif(strcmp(initialValues,"randn"))
        upLim = 1;
        loLim = -1;
        a = loLim + (upLim-loLim) .* randn(p, N);
    end
    grad_J = zeros(p,N);

    a_Steep_Conv = zeros(2,cycNumMax);

    c=0;
    cc(nn,1)=0;
    for nn = 1:N
        R(:,:,nn) = toeplitz(r_0(:,nn));
         
        
        % adding new method for steepest descent, I think it's wrong
         
         
         [~,lambda] = eig(R(:,:,nn));
         max_mu = 2/(max(diag(lambda)));
         mu = tuningMu * 2/(max(diag(lambda)));
         
         rapporto = max(diag(lambda)) / min(diag(lambda));
         
         % for ss = 1:numberOfCycles
         if (convergenceTest == 1)
            if (nn == 600)
                 a_Steep_Conv(1,1) = a(20,nn);
                 a_Steep_Conv(2,1) = a(21,nn);
            end
         end
         
         grad_J(:,nn) = -(2 *r_1(:,nn)) + (2* R(:,:,nn) * a(:,nn)); 
         a(:,nn) = a(:,nn) - (0.5 * mu *grad_J(:,nn));

         
        ccc = 2;
        
        while norm(grad_J(:,nn),2)> minThresh && ccc < cycNumMax % sqrt(p*0.0001)  
            if (convergenceTest == 1)
                if (nn == chosen_conv)
                     a_Steep_Conv(1,ccc) = a(comparing_conv_1,nn);
                     a_Steep_Conv(2,ccc) = a(comparing_conv_2,nn);
                end
            end
            R_a = (2* R(:,:,nn) * a(:,nn)); 
            grad_J(:,nn) = -(2 *r_1(:,nn)) + R_a;

            a(:,nn) = a(:,nn) - (0.5 * mu * ( grad_J(:,nn)));
            
            ccc=ccc+1;
        end
        if (convergenceTest == 1)
            if (nn == chosen_conv)
                a_Steep_Conv = a_Steep_Conv(:,1:ccc-1);
            end
        end
        if(strcmp(initialValues,"prev"))
            if(nn<N)
                a(:,nn+1) = a(:,nn);
            end
        end
        if verbose == 1 
            disp(nn)
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
    
    if (convergenceTest == 1) 
        a_test_Conv = lpc(st_signal(:,chosen_conv),p);
        figure 
        %for ii= 1:length(a_Steep_Conv)
        %    plot(-a_Steep_Conv(1,ii),-a_Steep_Conv(2,ii),'k-o')
        %    hold on
        %end
        plot(-a_Steep_Conv(1,:),-a_Steep_Conv(2,:),'k-o')
        hold on
        plot(a_test_Conv(1,comparing_conv_1 +1),a_test_Conv(1,comparing_conv_2 +1),'ro')
    end

    
    shapingFilters = H;
    
    whiteningFilters = A;
    
elseif strcmp(solveMode,"linSolve") 
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
    if verbose == 1 
        a_test = lpc(st_signal(:,ceil(N/3)),p);
        a_test = a_test';
        a_check = a(:,ceil(N/3));
        a_1_check = vertcat(1, -a_check);
        
        %[H_test,w] = freqz(1,a_1_test,"whole",M);
        figure('Renderer', 'painters', 'Position', [10 10 1000 600]);
        titlesize = 22;
        plot(a_test,"k+");
        hold on
        plot(a_1_check,"ro");
        legend('lpc function','custom function',Interpreter='latex');
        
        title("Test between LPC and custom made method",Interpreter='latex',FontSize=titlesize);
        
        A_fft = fft(a_1_check,M);
        H_fft = 1./A_fft;
        H_fft = H_fft';
        [H_freqz_lpc,w] = freqz(1,a_test,"whole",M);
        figure
        plot(w,abs(H_fft).^2,"k-*");
        hold on
        plot(w,abs(H(:,ceil(N/3))).^2,"b--");
        hold on
        plot(w,abs(H_freqz_lpc).^2,"ro");
        title("Test between LPC and custom made filters ");
        xlim([0 0.5]);
    end
    shapingFilters = H;
    
    whiteningFilters = A;
end