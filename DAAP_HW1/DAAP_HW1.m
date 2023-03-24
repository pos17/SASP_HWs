%==========================================================================
%                           DAAP HW1
%                           
%==========================================================================

% xcorr, circshift, cell

instr_name = "piano.wav";
instr_name_recon = "piano_recon.wav";
speech_name = "speech.wav";




[instr_t,instr_Fs] = audioread(instr_name);
[speech_t,speech_Fs] = audioread(speech_name);

instr_t_len =length(instr_t);

taps = 1000;
p = taps;
wl =  1024;  % windowLenght =

t_buckets = ceil(instr_t_len/wl);
instr_t(end+1:wl*t_buckets) = 0;
st_signal = reshape(instr_t,[wl t_buckets]);

r = zeros(p+1, t_buckets);

M = wl;
for nn = 1:t_buckets
    for ii_kk = 1:p+1
        for mm = 1:M-(ii_kk)
            r(ii_kk,nn) = r(ii_kk,nn)+ st_signal(mm)*st_signal(mm+ii_kk);
        end
    end 
end
r_0 = r(1:end-1,:);
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

%Hp = hamming(63,'periodic');
%wvtool(Hp)

%Hp = hamming(63,'symmetric');
%wvtool(Hs)