function withnoise = add_noise(speech,input,SNR,fs,fs_noise)
% speech is the clean signal without noise, and its sample frequency is fs;
% file is the noise's path and name, and the SNR is signal to noise ratio in dB. [input,fs_noise] = audioread(file);
noise_size = length(input);
speech_size = length(speech);
%compare different frequency of the speech and noise, if not equal,resample
if fs_noise~=fs
temp_noise = resample(input,fs,fs_noise);
else
temp_noise = input;
end
if noise_size >= speech_size
noise = temp_noise(1:speech_size);
else
noise = [temp_noise;zeros(speech_size-noise_size)];%if noise is short, add 0 at the end.
end
noise = noise-mean(noise);
signal_power = 1/speech_size*sum(speech.*speech); 
noise_variance = signal_power / (10^(SNR/10)); 
noise = sqrt(noise_variance)/std(noise)*noise; 
withnoise = speech + noise;
end