clear all
close all
clc

if not(isfolder("plots"))
    mkdir("Plots")
end

%% Import Input Audio Signal
[Vin,~] = audioread('ExpSweep.wav');

%% LTSpice Files for Ground-Truth
[OutLowSpice,~]=audioread('outlowsweep.wav');
[OutMidSpice,~]=audioread('outmidsweep.wav');
[OutHighSpice,FsLTSpice]=audioread('outhighsweep.wav');
TsLTSpice=1/FsLTSpice;

%% Sampling frequency (to be varied: FsLTSpice/downSampFact, with downSampFact={4,3,2})
downSampFact=1;
Fs = FsLTSpice/downSampFact; 

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
Z27 = RspkMid;
Z39 = RspkLow;

% Condenser adaptation

Z9 = Ts/(2*C1);
Z18 = Ts/(2*C2);
Z21 = Ts/(2*C3);
Z30 = Ts/(2*C4);
Z36 = Ts/(2*C5);
Z41 = Ts/(2*C6);

% Inductors adaptation

Z12 = (2*L1)/Ts;
Z13 = (2*L2)/Ts;
Z24 = (2*L3)/Ts;
Z31 = (2*L4)/Ts;


% Series and Parallels nodes 
% HIGH SECTION 

Z11 = (Z10*Z12) /(Z10+Z12);
Z7 = Z11;
Z8 = Z7 +Z9;


% MID SECTION

Z29 = Z28 +Z30;
Z26 = Z29;
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

Z1 = Z8;
Z4 = Z32;
Z6 = Z14;

Z5 = (Z4*Z6) /(Z4+Z6);
Z2 = Z5;
Z3 = (Z1*Z2) /(Z1+Z2);

%% Computation of Scattering Matrices 

%Sp1 = ParallelAdaptor(Z1,Z2,Z3);
% Matrix calculation
G1 = Z1^(-1);
G2 = Z2^(-1);
G3 = Z3^(-1);
portMatrix = zeros(3);
totG = sum([G1 G2 G3]);
onesCol = ones(3,1);
onesDiag = diag(onesCol);
rowG = [G1,G2,G3];
portMatrix = ((2/totG)*onesCol*rowG)-onesDiag;
Sp1 =portMatrix;

%Sp2 = ParallelAdaptor(Z4,Z5,Z6);
G1 = Z4^(-1);
G2 = Z5^(-1);
G3 = Z6^(-1);
portMatrix = zeros(3);
totG = sum([G1 G2 G3]);
onesCol = ones(3,1);
onesDiag = diag(onesCol);
rowG = [G1,G2,G3];
portMatrix = ((2/totG)*onesCol*rowG)-onesDiag;
Sp2 =portMatrix;

% HIGH

%Ss1 = SeriesAdaptor(Z7,Z8,Z9);
Zs1 = Z7;
Zs2 = Z8;
Zs3 = Z9;
totZ = sum([Zs1 Zs2 Zs3]);
onesRow = ones(1,3);
onesDiag = diag(onesRow);
colZ = [Zs1;Zs2;Zs3];
portMatrix = onesDiag - ((2/totZ)*(colZ)*onesRow);
Ss1 = portMatrix;


%Sp3 = ParallelAdaptor(Z10,Z11,Z12);
G1 = Z10^(-1);
G2 = Z11^(-1);
G3 = Z12^(-1);
totG = sum([G1 G2 G3]);
onesCol = ones(3,1);
onesDiag = diag(onesCol);
rowG = [G1,G2,G3];
portMatrix = ((2/totG)*onesCol*rowG)-onesDiag;
Sp3 =portMatrix;

% MID

%Ss2 = SeriesAdaptor(Z13,Z14,Z15);
Zs1 = Z13;
Zs2 = Z14;
Zs3 = Z15;
totZ = sum([Zs1 Zs2 Zs3]);
onesRow = ones(1,3);
onesDiag = diag(onesRow);
colZ = [Zs1;Zs2;Zs3];
portMatrix = onesDiag - ((2/totZ)*(colZ)*onesRow);
Ss2 = portMatrix;


%Sp4 = ParallelAdaptor(Z16,Z17,Z18);
G1 = Z16^(-1);
G2 = Z17^(-1);
G3 = Z18^(-1);
portMatrix = zeros(3);
totG = sum([G1 G2 G3]);
onesCol = ones(3,1);
onesDiag = diag(onesCol);
rowG = [G1,G2,G3];
portMatrix = ((2/totG)*onesCol*rowG)-onesDiag;
Sp4=portMatrix;

%Ss3 = SeriesAdaptor(Z19,Z20,Z21);
Zs1 = Z19;
Zs2 = Z20;
Zs3 = Z21;
totZ = sum([Zs1 Zs2 Zs3]);
onesRow = ones(1,3);
onesDiag = diag(onesRow);
colZ = [Zs1;Zs2;Zs3];
portMatrix = onesDiag - ((2/totZ)*(colZ)*onesRow);
Ss3 = portMatrix;

%Sp5 = ParallelAdaptor(Z22,Z23,Z24);
G1 = Z22^(-1);
G2 = Z23^(-1);
G3 = Z24^(-1);
portMatrix = zeros(3);
totG = sum([G1 G2 G3]);
onesCol = ones(3,1);
onesDiag = diag(onesCol);
rowG = [G1,G2,G3];
portMatrix = ((2/totG)*onesCol*rowG)-onesDiag;
Sp5 =portMatrix;

%Sp6 = ParallelAdaptor(Z25,Z26,Z27);
G1 = Z25^(-1);
G2 = Z26^(-1);
G3 = Z27^(-1);
portMatrix = zeros(3);
totG = sum([G1 G2 G3]);
onesCol = ones(3,1);
onesDiag = diag(onesCol);
rowG = [G1,G2,G3];
portMatrix = ((2/totG)*onesCol*rowG)-onesDiag;
Sp6 =portMatrix;

%Ss4 = SeriesAdaptor(Z28,Z29,Z30);
Zs1 = Z28;
Zs2 = Z29;
Zs3 = Z30;
totZ = sum([Zs1 Zs2 Zs3]);
onesRow = ones(1,3);
onesDiag = diag(onesRow);
colZ = [Zs1;Zs2;Zs3];
portMatrix = onesDiag - ((2/totZ)*(colZ)*onesRow);
Ss4 = portMatrix;

% LOW 

%Ss5 = SeriesAdaptor(Z31,Z32,Z33);
Zs1 = Z31;
Zs2 = Z32;
Zs3 = Z33;
totZ = sum([Zs1 Zs2 Zs3]);
onesRow = ones(1,3);
onesDiag = diag(onesRow);
colZ = [Zs1;Zs2;Zs3];
portMatrix = onesDiag - ((2/totZ)*(colZ)*onesRow);
Ss5 = portMatrix;

%Sp7 = ParallelAdaptor(Z34,Z35,Z36);
G1 = Z34^(-1);
G2 = Z35^(-1);
G3 = Z36^(-1);
portMatrix = zeros(3);
totG = sum([G1 G2 G3]);
onesCol = ones(3,1);
onesDiag = diag(onesCol);
rowG = [G1,G2,G3];
portMatrix = ((2/totG)*onesCol*rowG)-onesDiag;
Sp7 =portMatrix;

%Sp8 = ParallelAdaptor(Z37,Z38,Z39);
G1 = Z37^(-1);
G2 = Z38^(-1);
G3 = Z39^(-1);
portMatrix = zeros(3);
totG = sum([G1 G2 G3]);
onesCol = ones(3,1);
onesDiag = diag(onesCol);
rowG = [G1,G2,G3];
portMatrix = ((2/totG)*onesCol*rowG)-onesDiag;
Sp8 =portMatrix;

Ss6 = SeriesAdaptor(Z40,Z41,Z42);
Zs1 = Z40;
Zs2 = Z41;
Zs3 = Z42;
totZ = sum([Zs1 Zs2 Zs3]);
onesRow = ones(1,3);
onesDiag = diag(onesRow);
colZ = [Zs1;Zs2;Zs3];
portMatrix = onesDiag - ((2/totZ)*(colZ)*onesRow);
Ss6 = portMatrix;

%% Initialization of Waves
a = zeros(42,Nsamp+1);
b = zeros(42,Nsamp+1);

%% Initialize Output Signals
% Low
VoutLow=zeros(length(Vin),1);
% Mid
VoutMid=zeros(length(Vin),1);
% High
VoutHigh=zeros(length(Vin),1);

ii=0;
ii = ii+1; % to take into account matlab shift in elements 
while (ii<=Nsamp)
    ii=ii+1;

    %% Manage Dynamic Elements
    
    % HIGH 
    a(9,ii) = b(9,ii-1);
    a(10,ii) = 0;
    a(12,ii) = -b(12,ii-1);
    
    
    % MID 
    a(13,ii) = -b(13,ii-1);
    a(18,ii) = b(18,ii-1);
    a(21,ii) = b(21,ii-1);
    a(24,ii) = -b(24,ii-1);
    a(27,ii) = 0;
    a(28,ii) = 0;
    a(30,ii) = b(30,ii-1);
    
    %LOW
    a(31,ii) = -b(31,ii-1);
    a(36,ii) = b(36,ii-1);
    a(39,ii) = 0;
    a(41,ii) = b(41,ii-1);
    a(42,ii) = 0; 

    %% Forward Scan
    
    % HIGH 
    b(11,ii) = Sp3(2,:)*a(10:12,ii);
    a(7,ii) = b(11,ii);
    b(8,ii) = Ss1(2,:)*a(7:9,ii);
    
    % MID
    
    b(29,ii) = Ss4(2,:)*a(28:30,ii);
    a(26,ii) = b(29,ii);
    b(25,ii) = Sp6(1,:)*a(25:27,ii);
    a(23,ii) = b(25,ii);
    b(22,ii) = Sp5(1,:)*a(22:24,ii);
    a(20,ii) = b(22,ii);
    b(19,ii) = Ss3(1,:)*a(19:21,ii);
    a(17,ii) = b(19,ii);
    b(16,ii) = Sp4(1,:)*a(16:18,ii);
    a(15,ii) = b(16,ii);
    b(14,ii) = Ss2(2,:)*a(13:15,ii);
    
    %LOW 

    b(40,ii) = Ss6(1,:)*a(40:42,ii);
    a(38,ii) = b(40,ii);
    b(37,ii) = Sp8(1,:)*a(37:39,ii);
    a(35,ii) = b(37,ii);
    b(34,ii) = Sp7(1,:)*a(34:36,ii);
    a(33,ii) = b(34,ii);
    b(32,ii) = Ss5(2,:)*a(31:33,ii);
    
    % Main circuit parallels

   
    
    a(1,ii) = b(8,ii);
    a(4,ii) = b(32,ii);
    a(6,ii) = b(14,ii);
    b(5,ii) = Sp2(2,:)*a(4:6,ii);
    a(2,ii) = b(5,ii);
    b(3,ii) = Sp1(3,:)*a(1:3,ii);

    %% Local Root Scattering
    jj = ii-1; 

    a(3,ii) =  2* Vin(jj,1) - b(3,ii);

    %% Backward Scan
    b(1,ii) = Sp1(1,:)*a(1:3,ii);
    b(2,ii) = Sp1(2,:)*a(1:3,ii);
    a(8,ii) = b(1,ii);
    a(5,ii) = b(2,ii);
    b(4,ii) = Sp2(1,:)*a(4:6,ii);
    b(6,ii) = Sp2(3,:)*a(4:6,ii);
    a(14,ii) = b(6,ii);
    a(32,ii) = b(4,ii);

    % HIGH 
    
    b(7,ii) = Ss1(1,:)*a(7:9,ii);
    b(9,ii) = Ss1(3,:)*a(7:9,ii);
    a(11,ii) = b(7,ii);
    b(10,ii) = Sp3(1,:)*a(10:12,ii);
    b(12,ii) = Sp3(3,:)*a(10:12,ii);
    
    % MID

    b(13,ii) = Ss2(1,:)*a(13:15,ii);
    b(15,ii) = Ss2(3,:)*a(13:15,ii);
    a(16,ii) = b(15,ii);
    b(17,ii) = Sp4(2,:)*a(16:18,ii);
    b(18,ii) = Sp4(3,:)*a(16:18,ii);
    a(19,ii) = b(17,ii);
    b(20,ii) = Ss3(2,:)*a(19:21,ii);
    b(21,ii) = Ss3(3,:)*a(19:21,ii);
    a(22,ii) = b(20,ii);
    b(23,ii) = Sp5(2,:)*a(22:24,ii);
    b(24,ii) = Sp5(3,:)*a(22:24,ii);
    a(25,ii) = b(23,ii);
    b(26,ii) = Sp6(2,:)*a(25:27,ii);
    b(27,ii) = Sp6(3,:)*a(25:27,ii);
    a(29,ii) = b(26,ii);
    b(28,ii) = Ss4(1,:)*a(28:30,ii);
    b(30,ii) = Ss4(3,:)*a(28:30,ii);
    
    % LOW

    b(31,ii) = Ss5(1,:)*a(31:33,ii);
    b(33,ii) = Ss5(3,:)*a(31:33,ii);
    a(34,ii) = b(33,ii);
    b(35,ii) = Sp7(2,:)*a(34:36,ii);
    b(36,ii) = Sp7(3,:)*a(34:36,ii);
    a(37,ii) = b(35,ii);
    b(38,ii) = Sp8(2,:)*a(37:39,ii);
    b(39,ii) = Sp8(3,:)*a(37:39,ii);
    a(40,ii) = b(38,ii);
    b(41,ii) = Ss6(2,:)*a(40:42,ii);
    b(42,ii) = Ss6(3,:)*a(40:42,ii);
    
    
    %% Read Output
    
    
    % HIGH 

    VoutHigh(jj,1) = (a(10,ii) - b(10,ii))/2;
    

    % MID 

    VoutMid(jj,1) = (a(27,ii) + b(27,ii))/2;
    

    % LOW 

    VoutLow(jj,1) = (a(39,ii) - b(39,ii))/2;

    if(ii==floor(Nsamp/2)) 
        disp("debug")
    end


end

%%
V_3 = (a(3,:)+b(3,:))/2;
V_7 = (a(7,:)+b(7,:))/2;
V_8 = (a(8,:)+b(8,:))/2;
V_9 = (a(8,:)+b(9,:))/2;

figure
plot(V_3)
hold on 
plot(Vin)
figure
plot(V_7);
hold on 
plot(V_8);
hold on 
plot(V_9);
legend("V_7","V_8","V_9")

I_7 = (a(7,:)-b(7,:))/(2*Z7);
I_8 = (a(8,:)-b(8,:))/(2*Z8);
I_9 = (a(9,:)-b(9,:))/(2*Z9);


figure
plot(I_7);
hold on 
plot(I_8);
hold on 
plot(I_9);


%%

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
saveas(gcf,strcat("Plots/","Values Comparison - downsampling factor = " + num2str(downSampFact)),"png")

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
saveas(gcf,strcat("Plots/","Error plot - downsampling factor = " + num2str(downSampFact)),"png")

rms(VoutLow)