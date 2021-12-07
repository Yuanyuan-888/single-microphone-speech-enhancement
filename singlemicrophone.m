clc;
clear all;
[s fs] = audioread('clean_speech.wav');
[babble fbn] = audioread('babble_noise.wav');
[shaped fsn] = audioread('Speech_shaped_noise.wav');
[nonstat fan] = audioread('aritificial_nonstat_noise.wav');
snr_input =5;%10 dB
y = add_noise(s, nonstat, snr_input, fs,fan);

%sound(y,fs);% noisy signal

%% Framing & Apply transform
frame_length = 20*1e-3 * fs; % 20ms / 320samples
hop_length = 0.50 * frame_length; % 50% 
n = (length(y)-hop_length)/(frame_length-hop_length);
num_frame = floor(n)+1;% calculate the number of frames
y_t = zeros(frame_length, num_frame-1);
for idx = 1:(num_frame-1)
    front = 1+(idx-1)*(frame_length-hop_length);    
    rear = front + frame_length-1;
    y_t(:, idx) = y(front: rear);
end
last_segment = y(1+(num_frame-1)*(frame_length-hop_length):end);% consider the last segment
last_segment = [last_segment; zeros(frame_length-length(last_segment), 1)];
y_t = [y_t last_segment];
window = hann(frame_length);%apply hanning window
y_k = zeros(size(y_t));
for idx = 1 : num_frame
    y_k(:, idx) = fft(window .* y_t(:, idx));%FFT transform
end


%% Noisy speech PSD
L = 12; % number of periodogram segement 
Pyy = Bartlett(y_k,L);% Barlett noisy speech PSD estimate 

%% Noise PSD estimator: choose MS or MMSE
%% MS
M = 12; % num of segements
B = 1; % bias compensation
Pnn = NoisePSDMS(Pyy,M,B);
%% MMSE
%Pnn = NoisePSDMMSE(Pyy);


%% Target speech PSD estimator: choose ML or DD
%% ML 
SNR_ml = SpeechPSDML(Pyy,Pnn);
%% DD
%alpha = 0.98;%alpha 0.96-0.99
%SNR_dd = SpeechPSDDD(Pyy,Pnn,alpha);

%% Target Estimate:
%% power spectral substraction
%s_hat_k = spectral_substraction(Pyy, Pnn, y_k, 0.2);
%% Wiener gain 
s_hat_k = wiener(y_k,SNR_ml);

%% Inverse transform & Overlapp-add
s_t = ifft(s_hat_k);% Inverse FFT transform
s_t_est1 = s_t(1:frame_length-hop_length/2, 1);
s_t_est2 = s_t(1+hop_length/2:end,end);
s_t(:, 1) = [];
s_t(:, end) = [];
s_t(1:hop_length/2, :) = [];
s_t(end-(hop_length/2-1): end, :) = [];
s_t = reshape(s_t, [], 1);
s_t_est = [s_t_est1; s_t;s_t_est2];
s_t = real(s_t);


%% Evaluation
length=min(size(s),size(s_t));
%% SNR
% n = s_t(1:length)-s(1:length);%after noise reduction, noise
% P_clean=norm(s(1:length)).^2;
% P_yynow=norm(s_t(1:length)).^2;
% Pnn_after_method1 = P_yynow-P_clean;
% %P_clean = sum(abs(s(1:length)).^2);
% Pnn_after=sum(abs(n(1:length)).^2);
% snr_after = 10.*log10(P_clean/Pnn_after);
% SNR_improve=snr_after-snr_input;
% %SNR_improve=SNR_ml-snr_input;
frame_len = 0.02 * fs; %20ms 
overlap_len = 0.5 * frame_len;
Sl = ifft(s_hat_k);
Sl = Sl(1:frame_len,:);
Sn = overlap_add(Sl,frame_len,overlap_len); 
Sn = real(Sn(1:length));
%% listening test
%soundsc(clean,fs);
%soundsc(Sn,fs);
%% segmental SNR
overlap_len = 0;
frame_num = floor((length-overlap_len)/(frame_len-overlap_len))+1; 
for i = 1:(frame_num-1)
start = (i-1)*(frame_len-overlap_len) + 1; 
stop = start + frame_len - 1;
y_c(:, i) = s(start:stop);
end
last_original = s((frame_num-1)*(frame_len-overlap_len) + 1:length); 
last = [last_original; 
zeros(frame_len-length(last_original),1)];
y_c = [y_c last];
overlap_len = 0;
for i = 1:(frame_num-1)
start = (i-1)*(frame_len-overlap_len) + 1; 
stop = start + frame_len - 1;
S(:, i) = Sn(start:stop);
end
last_original = Sn((frame_num-1)*(frame_len-overlap_len) + 1:length); 
last = [last_original; 
zeros(frame_len-length(last_original),1)];
S = [S last];
sum_Pyy = sum(abs(S).^2);
% Sc = clean;
nl = S-y_c;
sum_Pnn = sum(abs(nl).^2);
snr = 10.*log(sum_Pyy./sum_Pnn);

for i=1:1804
    
   if snr(:,i)>-10
        if snr(:,i)<35
            snr(:,i)=snr(:,i);
        else
                snr(:,i)=35;
        end
    end
end               
sum_snr = sum(snr);
snrseg = (1/1804).*sum_snr; 
snr_imprv = snrseg-snr_input;
%% STOI
d = stoi(s(1:length), s_t(1:length), fs);

figure;
subplot(321);
plot(s(10000:170000));
xlim([10000 170000]);
title('clean speech');
xlabel('Sample')
ylabel('Amplitude')

subplot(322);
spectrogram(s(10000:170000),window,hop_length,frame_len,fs,'yaxis');
colormap; 
colormap(jet);
title('Clean Speech');

subplot(323);
plot(y(10000:170000));
xlim([10000 170000]);
title('noisy speech');
xlabel('Sample')
ylabel('Amplitude')

subplot(324);
spectrogram(y(10000:170000),window,hop_length,frame_len,fs,'yaxis');
colormap; 
colormap(jet);
title('Noisy Speech');

subplot(325);
plot(s_t(10000:170000));
xlim([10000 170000]);
title('filtered speech ');
xlabel('Sample')
ylabel('Amplitude')


subplot(326); 
spectrogram(s_t(10000:170000),window,hop_length,frame_len,fs,'yaxis');
colormap; 
colormap(jet);
title('Enhanced Speech');

