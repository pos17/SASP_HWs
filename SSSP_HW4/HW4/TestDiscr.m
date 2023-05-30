
jw = linspace(0,pi,240000);
% Ts = 1/48000
% jww3 = 2/Ts * tan(jwtilde* Ts/2)
figure
%plot(jw,jw-jw)
%hold on

lim = 24000;
Fs = 192000;
Ts = 1/Fs;
jwtilde = linspace(0,pi,Fs/2+1);
jww = 2/Ts * tan(jwtilde* Ts/2);
f = jww/(2*pi)*Fs;
jw_minus = jwtilde; %/(2*pi)*Fs;
f_minus = jwtilde/(2*pi)*Fs;
%found = find(f_minus==lim);
f_minus = f_minus(1:Fs/8+1);
jwtilde =jwtilde(length(f_minus));
f=f(length(f_minus));
%plot(jwtilde,abs(jww -jw_minus))
plot(jwtilde,abs(f -f_minus))

hold on

Fs = 96000;
Ts = 1/Fs;
jwtilde = linspace(0,pi,Fs/2+1);
jww = 2/Ts * tan(jwtilde* Ts/2);
f = jww/(2*pi)*Fs;
jw_minus =  jwtilde; %/(2*pi)*Fs;
f_minus = jwtilde/(2*pi)*Fs;
f_minus = f_minus(1:Fs/4+1);
jwtilde =jwtilde(length(f_minus));
f=f(length(f_minus));
%plot(jwtilde,abs(jww -jw_minus))
plot(jwtilde,abs(f -f_minus))
hold on

Fs = 48000;
Ts = 1/Fs;
jwtilde = linspace(0,pi,Fs/2+1);
jww = 2/Ts * tan(jwtilde* Ts/2);
f = jww/(2*pi)*Fs;
jw_minus =  jwtilde; %/(2*pi)*Fs
f_minus = jwtilde/(2*pi)*Fs;
f_minus = f_minus(1:Fs/2+1);
jwtilde =jwtilde(length(f_minus));
f=f(length(f_minus));
%plot(jwtilde,abs(jww -jw_minus))
plot(jwtilde,abs(f -f_minus))
grid minor
%ylim([-1, max(jww3-jwtilde+1)])

legend("Continuous","Fs=192k","Fs=96k","Fs=48k")