function psd_b = Bartlett (y, M)
%y is the signal
%M is the number of periodogram segement to average
%output is the power spectrum density estimated by Bartlett method
[L num] = size(y); % length of segements & num of segements
periodogram = zeros(L, num);% periodogram of each segement
for idx = 1:num
    periodogram(:,idx) = abs(y(:,idx)) .^2;
end
psd_b = zeros(L, num);
for l = 1: num% Bartlett PSD is the average of M periodogram segments
    if l<M
        psd_b(:,l) = 1/l*sum(periodogram(:, 1:l), 2);
    else
        psd_b(:,l) = 1/M * sum(periodogram(:, l-M+1:l), 2);
    end
end