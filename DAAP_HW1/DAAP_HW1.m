%==========================================================================
%                           DAAP HW1
%                           
%==========================================================================

% xcorr, circshift, cell

instr_name = "piano.wav";
instr_name_recon = "piano_recon.wav";
speech_name = "speech.wav";
st_res_name = "whitening_filter_res.wav";



[instr_t,instr_Fs] = audioread(instr_name);
[speech_t,speech_Fs] = audioread(speech_name);

instr_t_len =length(instr_t);

taps = 220;
p = taps;
wl =  1024;  % windowLenght =

t_buckets = ceil(instr_t_len/wl);
instr_t(end+1:wl*t_buckets) = 0;
st_signal = reshape(instr_t,[wl t_buckets]);

r = zeros(p+1, t_buckets);

M = wl;
for nn = 1:t_buckets
    for ii_kk = 0:p
        for mm = 1:M-(ii_kk)
            r(ii_kk+1,nn) = r(ii_kk+1,nn)+ st_signal(mm)*st_signal(mm+ii_kk);
        end
    end 
end
r_0 = r(1:end-1,:);
plot(r_0);
r_1 = r(2:end,:);
R = zeros(p,p,t_buckets);
a = zeros(p, t_buckets);
for nn = 1:t_buckets
    R(:,:,nn) = toeplitz(r_0(:,nn));
    a(:,nn) = R(:,:,nn)\r_1(:,nn);
end


R_1 = R(:,:,1);
%% Reconstruction of the prediction using prediction filters 

reconstructed_Sig = zeros(wl, t_buckets);

for nn = 1:t_buckets
    for mm = p+1:M
        for kk = 1:p
            reconstructed_Sig(mm,nn) = reconstructed_Sig(mm,nn) ...
                + a(kk,nn)*st_signal(mm-kk+1,nn);
        end
    end
end
reconstructed_Sig = reshape(reconstructed_Sig,[t_buckets *wl 1]);

audiowrite(instr_name_recon,reconstructed_Sig,instr_Fs);

%% Z_transform of the filter represented as spectrum of the filter
n = wl;

w = linspace(-Fs/2,Fs/2,n);

A = zeros(length(w),t_buckets);
z= exp(1i*w);
sum_A =zeros(length(w),t_buckets);
for nn=1:t_buckets
    for kk = 1:p
        sum_A(:,nn) = a(kk,nn)*z.^(-kk);
    end
end
A = 1- sum_A;
Aishift = zeros(wl,t_buckets);
for nn = 1: t_buckets
    Aishift(:,nn) = ifftshift(A(:,nn));
end 
A1 = Aishift(1:512,:);
%% STFT of the signal 
w1   = linspace(0,instr_Fs/2,wl/2);
st_signal_w = zeros(wl,t_buckets);
for nn = 1:t_buckets
    st_signal_w(:,nn) = fft(st_signal(:,nn));
end

Fs = instr_Fs;

t =linspace(0,instr_t_len/instr_Fs,t_buckets);
st_signal_w1 = st_signal_w(1:512,:);
surf(t,w1,abs(st_signal_w1),EdgeColor="none");
ylim([0 5000])
view(0,90)

%% signal filtering



st_res =  zeros(wl, t_buckets); 
st_res_w = zeros(wl,nn);
for nn = 1:t_buckets
    st_res_w(:,nn) = A(:,nn) .* fftshift(st_signal_w(:,nn));
    st_res_w(:,nn) = ifftshift(st_res_w(:,nn));
    st_res(:,nn) = ifft(st_res_w(:,nn));
end

st_res_lin = reshape(st_res,[t_buckets *wl 1]);
audiowrite(st_res_name,abs(st_res_lin),instr_Fs);


%% plot of a random element of the matrix for the signal and for the filter 
plot(w1,abs(st_signal_w1(:,400)))
hold on 
plot(w1,abs(A1(:,400)))
xlim([0,500])

%Hp = hamming(63,'periodic');
%wvtool(Hp)

%Hp = hamming(63,'symmetric');
%wvtool(Hs)