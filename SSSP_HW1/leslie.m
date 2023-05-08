function [y,y_lpf,y_hpf,y_hp_sdf] = leslie(x, Fs, freq)
%Leslie Speaker Emulation
%
% J. Pekonen et al. Computationally Efficient Hammond Organ Synthesis
% in Proc. of the 14th International Conference on Digital Audio
% Effects(DAFx-11), Paris, France, Sept. 19-23, 2011

% length of the input signal
N = length(x);

% global modulator parameters
alpha=0.9;
% tremble spectral delay filter parameter 
Ms_t=0.2;
Mb_t=-0.75;
N_sdf_t=4;
% bass spectral delay filter parameter 
Ms_b=0.04;
Mb_b=-0.92;
N_sdf_b=3;

% cross-over network design
fc=800;                 % cutoff frequency

%TODO: compute the coefficients for the two 4th order butterworth filters
%with cutoff frequency fc
N_bw = 4;
[b_lp, a_lp] = butter(N_bw,fc/(Fs/2)); %LPF design
[b_hp, a_hp] = butter(N_bw,fc/(Fs/2),'high'); %HPF design

Hd_hp = dfilt.df1(b_hp,a_hp)

Hd_lp = dfilt.df1(b_lp,a_lp)
a_lp1 = a_lp(2:end)';
a_hp1 = a_hp(2:end)';
b_lp1 = b_lp';
b_hp1 = b_hp';


% crossover plot
figure()
freqz(b_lp,a_lp,[],Fs)
subplot(2,1,1)
ylim([-100 20])
xlim([0 2000])
legend("lpf")
grid on 
figure()
freqz(b_hp,a_hp,[],Fs)
ylim([-100 20])
xlim([0 2000])
grid on 
legend("hpf")

% allocate input and output buffers for IIR filters
% hp filter buffers
hpf.state=zeros(N_bw,1);
hpf.in=zeros(N_bw+1,1);
% lp filter buffers
lpf.state=zeros(N_bw,1);
lpf.in=zeros(N_bw+1,1);
% treble sdf filter buffers
sdf_h.state=zeros(N_sdf_t+1,1);
sdf_h.in=zeros(N_sdf_t+1,1);
% bass sdf filter buffers
sdf_b.state=zeros(N_sdf_b+1,1);
sdf_b.in=zeros(N_sdf_b+1,1);

% modulators
 m_b = 0;%$Ms_b*sin(freq*n)+Mb_b % bass modulator
 m_t = 0;%Ms_t*sin((freq+0.1)*n)+Mb_t % tremble modulator


% for n=1:2000
%     m_t(1,n) = Ms_t*sin((freq+0.1)*n)+Mb_t;
% end 
% figure()
% n=1:1:2000;
% plot(n,m_t);
% grid on

%sample processing

y_hpf = zeros(N,1);
y_lpf = zeros(N,1);
y = zeros(N,1);
y_lp_sdf= zeros(N,1);
y_hp_sdf = zeros(N,1);
%lpf.in(1,1) = x(1,1);
%hpf.in(1,1) = x(1,1);
for n=1:N
    
    %index calculation for circular array
    nn = n-1;

    lpf.in(mod(nn,N_bw+1)+1,1) = x(n,1);
    % compute crossover network filters outputs
    for ii = 0:N_bw-1 
        y_lpf(n,1)= y_lpf(n,1) - (a_lp1(ii+1,1)*lpf.state(mod(nn-ii-1,N_bw)+1,1));
    end
    for ii = 0:N_bw 
        y_lpf(n,1)= y_lpf(n,1) + (b_lp1(ii+1,1)*lpf.in(mod(nn-ii,N_bw+1)+1,1));
    end
    lpf.state(mod(nn,N_bw)+1,1) = y_lpf(n,1);
    
    hpf.in(mod(nn,N_bw+1)+1,1) = x(n,1);
    % compute crossover network filters outputs
    for ii = 0:N_bw-1 
        y_hpf(n,1)= y_hpf(n,1) - (a_hp1(ii+1,1)*hpf.state(mod(nn-ii-1,N_bw)+1,1));
    end
    for ii = 0:N_bw 
        y_hpf(n,1)= y_hpf(n,1) + (b_hp1(ii+1,1)*hpf.in(mod(nn-ii,N_bw+1)+1,1));
    end
    hpf.state(mod(nn,N_bw)+1,1) = y_hpf(n,1);
    
    %if( y_hpf(n,1) ~= filter(b_hp,a_hp,x))

    %y_hpf=...

    % compute bass SDF output
    
    m_b = Ms_b*sin((freq/(Fs))*2*pi*n)+Mb_b;
    
    sdf_b.in(mod(nn,N_sdf_b+1)+1,1) = y_lpf(n,1);
    sdf_b.state(mod(nn,N_sdf_b+1)+1,1) = 0;
    for ii = 0:N_sdf_b 
        y_lp_sdf(n)= y_lp_sdf(n) +(nchoosek(N_sdf_b,ii)* (m_b^ii) * (sdf_b.in(mod(nn-(N_sdf_b-ii),N_sdf_b+1)+1,1) - sdf_b.state(mod(nn-ii,N_sdf_b+1)+1,1))); % - a_lp(ii+1)*lpf.state(mod(nn-ii-1,4)+1) + b_lp(ii+1)*lpf.in(mod(nn-ii,4)+1);
    end
    sdf_b.state(mod(nn,N_sdf_b+1)+1,1) = y_lp_sdf(n);
    
    
    
    %y_lp_sdf=...

    % compute treble SDF output
    
    m_t = Ms_t*sin(((freq+0.1)/(Fs))*2*pi*n)+Mb_t; % tremble modulator
    
    sdf_h.in(mod(nn,N_sdf_t+1)+1,1) = y_hpf(n,1);
    sdf_h.state(mod(nn,N_sdf_t+1)+1,1) = 0;
    for ii = 0:N_sdf_t 
        y_hp_sdf(n)= y_hp_sdf(n) +(nchoosek(N_sdf_t,ii)* m_t^ii * (sdf_h.in(mod(nn-(N_sdf_t-ii),N_sdf_t+1)+1,1) - sdf_h.state(mod(nn-ii,N_sdf_t+1)+1,1))); % - a_lp(ii+1)*lpf.state(mod(nn-ii-1,4)+1) + b_lp(ii+1)*lpf.in(mod(nn-ii,4)+1);
    end
    sdf_h.state(mod(nn,N_sdf_t+1)+1,1) = y_hp_sdf(n);
    %y_hp_sdf=...

    % implement AM modulation blocks
    y_lp_am= (1+alpha*m_b)*y_lp_sdf(n);
    y_hp_am= (1+alpha*m_t)*y_hp_sdf(n);

    y(n,1)= y_lp_am + y_hp_am;
    
end
end

