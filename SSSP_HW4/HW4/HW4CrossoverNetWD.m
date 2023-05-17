clear all
close all
clc

%% Import Input Audio Signal
[Vin,~] = audioread('ExpSweep.wav');

%% LTSpice Files for Ground-Truth
[OutLowSpice,~]=audioread('outlowsweep.wav');
[OutMidSpice,~]=audioread('outmidsweep.wav');
[OutHighSpice,FsLTSpice]=audioread('outhighsweep.wav');
TsLTSpice=1/FsLTSpice;

%% Sampling frequency (to be varied: FsLTSpice/downSampFact, with downSampFact={4,3,2})
downSampFact=2;
Fs =FsLTSpice/downSampFact; 

%% Downsample Input Signal
Vin=Vin([1:downSampFact:end]);

%% Sampling Period
Ts=1/Fs;
%% Number of Samples
Nsamp=length(Vin);
%% Simulated time
tstop=Nsamp*Ts;
%% Parameters of Dynamic Element
L1=0.35*10^(-3);
L2=0.35*10^(-3);
L3=3.5*10^(-3);
L4=3.5*10^(-3);
C1= 2.8*10^(-6);
C2= 2.8*10^(-6);
C3= 28*10^(-6);
C4= 4.7*10^(-6);
C5=28*10^(-6);
C6=47*10^(-6);
%% Resistive Parameters
R1=10;
RspkLow=8;
R2=10;
RspkMid=8;
RspkHigh=8;
%% WDF setting of free parameters (adaptation conditions)

% Resistences adaptation

Z28 = R1;
Z42 = R2;
Z10 = RspkHigh;
Z26 = RspkMid;
Z39 = RspkLow;

% Condenser adaptation

Z9 = Ts/(2*C1);
Z18 = Ts/(2*C2);
Z21 = Ts/(2*C3);
Z30 = Ts/(2*C4);
Z36 = Ts/(2*C5);
Z41 = Ts/(2*C6);

% Inductors adaptation

Z12 = 2*L1/Ts;
Z13 = 2*L2/Ts;
Z24 = 2*L3/Ts;
Z31 = 2*L4/Ts;


% Series and Parallels nodes 
% HIGH SECTION 

Z11 = (Z10*Z12) /(Z10+Z12);
Z7 = Z11;
Z8 = Z7 +Z9;


% MID SECTION

Z29 = Z28 +Z30;
Z27 = Z29;
Z25 = (Z26*Z27) /(Z26+Z27);
Z23 = Z25;
Z22 = (Z23*Z24) /(Z23+Z24);
Z20 = Z22;
Z19 = Z20 + Z21;
Z17 = Z19;
Z16 = (Z17*Z18) /(Z17+Z18);
Z15 = Z16;
Z14 = Z13 + Z15;


% LOW SECTION

Z40 = Z41 + Z42;
Z38 = Z40;
Z37 = (Z38*Z39) /(Z38+Z39);
Z35 = Z37;
Z34 = (Z35*Z36) /(Z35+Z36);
Z33 = Z34;
Z32 = Z31 + Z33;

% Main circuit parallels

Z1=Z8;
Z4 = Z32;
Z6 = Z14;

Z5 = (Z4*Z6) /(Z4+Z6);
Z2 = Z5;
Z3 = (Z1*Z2) /(Z1+Z2);

%% Computation of Scattering Matrices 

Sp1 = ParallelAdaptor(Z1,Z2,Z3);
Sp3 = ParallelAdaptor(Z4,Z5,Z6);

% HIGH

Ss1 = SeriesAdaptor(Z7,Z8,Z9);
Sp3 = ParallelAdaptor(Z10,Z11,Z12);

% MID

Ss2 = SeriesAdaptor(Z13,Z14,Z15);
Sp4 = ParallelAdaptor(Z16,Z17,Z18);
Ss3 = SeriesAdaptor(Z19,Z20,Z21);
Sp5 = ParallelAdaptor(Z22,Z23,Z24);
Sp6 = ParallelAdaptor(Z25,Z26,Z27);
Ss4 = SeriesAdaptor(Z28,Z29,Z30);

% LOW 

Ss5 = SeriesAdaptor(Z31,Z32,Z33);
Sp7 = ParallelAdaptor(Z34,Z35,Z36);
Sp8 = ParallelAdaptor(Z37,Z38,Z39);
Ss6 = SeriesAdaptor(Z40,Z41,Z42);

%% Initialization of Waves
a = zeros(42,2);
b = zeros(42,2);

%% Initialize Output Signals
% Low
VoutLow=zeros(size(Vin));
% Mid
VoutMid=zeros(size(Vin));
% High
VoutHigh=zeros(size(Vin));

ii=0;
while (ii<Nsamp)
    ii=ii+1;
    
    in = mod(ii,2) +1;
    in_1 = mod(ii-1,2) +1;

    %% Manage Dynamic Elements
    
    % HIGH 
    a(9,in) = b(9,in_1);
    a(10,in) = 0;
    a(12,in) = -b(12,in_1);
    
    % MID 
    a(13,in) = -b(13,in_1);
    a(18,in) = b(18,in_1);
    a(21,in) = b(21,in_1);
    a(24,in) = -b(24,in_1);
    a(26,in) = 0;
    a(28,in) = 0;
    a(30,in) = b(30,in_1);
    
    %LOW
    a(31,in) = -b(31,in_1);
    a(36,in) = b(36,in_1);
    a(39,in) = 0;
    a(41,in) = b(41,in_1);
    a(42,in) = 0;

    %% Forward Scan
    
    % HIGH 
    b(11,in) = Sp3(2,:)*a(10:12,in);
    a(7,in) = b(11,in);
    b(8,in) = Ss1(2,:)*a(7:9,in);
    
    % MID
    
    b(29,in) = Ss4(2,:)*a(28:30,in);
    a(26,in) = b(29,in);
    b(25,in) = Sp6(1,:)*a(25:27,in);
    a(23,in) = b(25,in);
    b(22,in) = Sp5(1,:)*a(22:24,in);
    a(20,in) = b(22,in);
    b(19,in) = Ss3(1,:)*a(19:21,in);
    a(17,in) = b(19,in);
    b(16,in) = Sp4(1,:)*a(16:18,in);
    a(15,in) = b(16,in);
    b(14,in) = Ss2(2,:)*a(13:15,in);
    
    %LOW 

    b(40,in) = Ss6(1,:)*a(40:42,in);
    a(38,in) = b(40,in);
    b(37,in) = Sp8(1,:)*a(37:39,in);
    a(35,in) = b(37,in);
    b(34,in) = Sp7(1,:)*a(34:36,in);
    a(33,in) = b(34,in);
    b(32,in) = Ss5(2,:)*a(31:33,in);
    
    % Main circuit parallels
    a(1,in) = b(8,in);
    a(4,in) = b(32,in);
    a(6,in) = b(14,in);
    a(2,in) = b(5,in);
    b(5,in) = Sp2(2,:)*a(4:6,in);
    b(3,in) = Sp1(3,:)*a(1:3,in);
    

    %% Local Root Scattering
    
    a(3,in) = 2* Vin(in,1) - b(3,in); 

    %% Backward Scan
    

    %% Read Output
    
    
end


%% Output Plots
figure
set(gcf, 'Color', 'w');
subplot(311)
plot(TsLTSpice*[1:length(OutLowSpice)],OutLowSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutLow,'b--','Linewidth',1); grid on; xlim([0,tstop]); 
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outLow}}$ [V]','Fontsize',16,'interpreter','latex');
legend('LTspice','WDF','Fontsize',16,'interpreter','latex');
title('Output Signals','Fontsize',18,'interpreter','latex');
subplot(312)
plot(TsLTSpice*[1:length(OutMidSpice)],OutMidSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutMid,'b--','Linewidth',1); grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outMid}}$ [V]','Fontsize',16,'interpreter','latex');
subplot(313)
plot(TsLTSpice*[1:length(OutHighSpice)],OutHighSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutHigh,'b--','Linewidth',1); grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outHigh}}$ [V]','Fontsize',16,'interpreter','latex');

%% Error Plots
figure
set(gcf, 'Color', 'w');
subplot(311)
plot(Ts*[1:Nsamp],OutLowSpice([1:downSampFact:end])-VoutLow,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outLow}}$ [V]','Fontsize',16,'interpreter','latex');
title(['Error Signals. $F_{\mathrm{s}}=$ ',num2str(Fs),' Hz'],'Fontsize',18,'interpreter','latex');
subplot(312)
plot(Ts*[1:Nsamp],OutMidSpice([1:downSampFact:end])-VoutMid,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outMid}}$ [V]','Fontsize',16,'interpreter','latex');
subplot(313)
plot(Ts*[1:Nsamp],OutHighSpice([1:downSampFact:end])-VoutHigh,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outHigh}}$ [V]','Fontsize',16,'interpreter','latex');

