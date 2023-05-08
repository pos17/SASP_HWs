function [y,y_lpf,y_hpf,y_hp_sdf] = leslie_o(x, Fs, freq)
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
[b_lp, a_lp] = butter(4,fc/(Fs/2)); %LPF design
[b_hp, a_hp] = butter(4,fc/(Fs/2),'high'); %HPF design

% crossover plot
figure()
freqz(b_lp,a_lp,[],Fs)
subplot(2,1,1)
ylim([-100 20])
xlim([0 2000])
hold on
grid on 
freqz(b_hp,a_hp,[],Fs)
ylim([-100 20])
xlim([0 2000])

% allocate input and output buffers for IIR filters
% hp filter buffers
buffer = 512;
hpf.state=zeros(buffer,1);
hpf.in=zeros(buffer,1);
% lp filter buffers
lpf.state=zeros(buffer,1);
lpf.in=zeros(buffer,1);
% treble sdf filter buffers
sdf_h.state=zeros(buffer,1);
sdf_h.in=zeros(buffer,1);
% bass sdf filter buffers
sdf_b.state=zeros(buffer,1);
sdf_b.in=zeros(buffer,1);

% modulators
% m_b(1,n) = Ms_b*sin(freq*n)+Mb_b % bass modulator
% m_t(1,n) = Ms_t*sin((freq+0.1)*n)+Mb_t % tremble modulator
% 

for n=1:2000
    m_t(1,n) = Ms_t*sin((freq+0.1)*n)+Mb_t;
end 
figure()
n=1:1:2000;
plot(n,m_t);
grid on

%sample processing
% for n=1:N
% 
%     % compute crossover network filters outputs
%     y_lpf=
% 
%     y_hpf=...
% 
%     % compute bass SDF output
%     y_lp_sdf=...
% 
%     % compute treble SDF output
%     y_hp_sdf=...
% 
%     % implement AM modulation blocks
%     y_lp_am=...
%     y_hp_am=...
% 
%     y(n)=...
% 
% end

end

