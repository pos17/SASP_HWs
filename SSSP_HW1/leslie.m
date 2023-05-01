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
[b_lp, a_lp]=... %LPF design
[b_hp, a_hp]=...  %HPF design

% allocate input and output buffers for IIR filters
% hp filter buffers
hpf.state=zeros(...,1);
hpf.in=zeros(...,1);
% lp filter buffers
lpf.state=zeros(...,1);
lpf.in=zeros(...,1);
% treble sdf filter buffers
sdf_h.state=zeros(...,1);
sdf_h.in=zeros(...,1);
% bass sdf filter buffers
sdf_b.state=zeros(...,1);
sdf_b.in=zeros(...,1);

% modulators
m_b=... % bass modulator
m_t=... % tremble modulator

%sample processing
for n=1:N

    % compute crossover network filters outputs
    y_lpf=...

    y_hpf=...

    % compute bass SDF output
    y_lp_sdf=...

    % compute treble SDF output
    y_hp_sdf=...

    % implement AM modulation blocks
    y_lp_am=...
    y_hp_am=...

    y(n)=...

end

end

