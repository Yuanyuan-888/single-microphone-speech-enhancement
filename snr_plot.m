function snr_impro = snr_plot(speech,noise,snr) 
[clean, fs] = audioread(speech);
snr_input = snr;
y = add_noise(clean, noise, snr_input, fs); frame_len = 0.02 * fs; %20ms
overlap_len = 0.5 * frame_len;
frame_num = floor((length(y)-overlap_len)/(frame_len-overlap_len))+1; y_l = zeros(frame_len,frame_num-1);
for i = 1:(frame_num-1)
start = (i-1)*(frame_len-overlap_len) + 1; stop = start + frame_len - 1;
y_l(:, i) = y(start:stop);
end
last_original = y((frame_num-1)*(frame_len-overlap_len) + 1:length(y)); last = [last_original; zeros(frame_len-length(last_original),1)];
y_l = [y_l last];
window = hann(frame_len);
for i = 1 : frame_num
y_k(:,i) = fft(window .* y_l(:,i),512);
end
%% use Bartlett estimate to compute Pyy
M = 12;
Pyy_p = zeros(size(y_k)); for i = 1 : frame_num
Pyy_p(:,i) = 1/frame_num * abs(y_k(:,i)) .^2;
end
Pyy = zeros(size(y_k));
for i = 1 : frame_num if i<M
Pyy(:,i) = 1/i * sum(Pyy_p(:,1:i),2);
else
Pyy(:,i) = 1/M * sum(Pyy_p(:,i-M+1:i),2);
    end
end
%% compute Pnn
%Pnn = min_sta(Pyy,M,1); Pnn = mmse(Pyy,20,0.8); %% esitimate SNR
SNR_ml = snr_ml(Pyy,Pnn); SNR_dd = snr_dd(Pyy,Pnn,0.97);
%% wiener smoother
Hk = SNR_ml./(SNR_ml+1);
Sk_hat = Hk .*abs(y_k).*exp(1j*angle(y_k));
Sl = ifft(Sk_hat);
Sl = Sl(1:frame_len,:);
Sn = overlap_add(Sl,frame_len,overlap_len); Sn = real(Sn(1:length(clean)));
%% segmental SNR
overlap_len = 0;
frame_num = floor((length(y)-overlap_len)/(frame_len-overlap_len))+1; for i = 1:(frame_num-1)
start = (i-1)*(frame_len-overlap_len) + 1; stop = start + frame_len - 1;
y_c(:, i) = clean(start:stop);
end
last_original = clean((frame_num-1)*(frame_len-overlap_len) + 1:length(y)); last = [last_original; zeros(frame_len-length(last_original),1)];
y_c = [y_c last];
overlap_len = 0;
for i = 1:(frame_num-1)
start = (i-1)*(frame_len-overlap_len) + 1; stop = start + frame_len - 1;
S(:, i) = Sn(start:stop);
end
last_original = Sn((frame_num-1)*(frame_len-overlap_len) + 1:length(y)); last = [last_original; zeros(frame_len-length(last_original),1)];
S = [S last];
sum_Pyy = sum(abs(S).^2);
nl = S-y_c;
sum_Pnn = sum(abs(nl).^2);
snr = 10.*log(sum_Pyy./sum_Pnn);
for i=1:1806
if snr(:,i)>-10
if snr(:,i)<35 snr(:,i)=snr(:,i);
else
snr(:,i)=35;
end
else snr(:,i)=-10;
end
end
sum_snr = sum(snr);
snrseg = 1/1086.*sum_snr; snr_impro = snrseg-snr_input; 
end