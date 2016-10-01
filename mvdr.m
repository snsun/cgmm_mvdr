% Created on 2016-08-19
% Author: Zhang Binbin
% About: MVDR matlab code
% Modified by Sining Sun 

[pcm, fs] = audioread('4.wav');
[num_point, num_channel] = size(pcm);
num_stat = 129; % 快拍数 % 2^M
tdoa_window = 400; % 估计TDOA window大小, 2^N
stft_len = 512;
frm_shift = 160;
frm_num = floor((num_point - tdoa_window) / frm_shift + 1);
output = zeros(num_point, 1);


spectrum = zeros(num_stat, num_channel, stft_len );
%space = round(tdoa_window / num_stat);
for s = 1:num_stat
    sub_data = pcm((s-1)*frm_shift+1:(s-1)*frm_shift+400, :).* repmat(hamming(tdoa_window), 1, num_channel);
    spectrum(s, :, :) = fft(sub_data, stft_len)';
end
corrvar = zeros(num_channel, num_channel, stft_len /2 + 1, num_stat);
for j = 1:num_stat
    for k = 1 : stft_len / 2 + 1
        % corrvar(:,:,k) = cov(spectrum(:, :, k));
        corrvar(:,:,k, j) = spectrum(j, :, k)' * conj(spectrum(j, :, k));
        corrvar(:, :, k, j) = corrvar(:, :, k, j) / trace(corrvar(:, :, k, j));
    end
end
    corrvar = mean(corrvar, 4);;
    
    
for j = 1:frm_shift:frm_shift*frm_num
    % decide tdoa window size
    segment_size = tdoa_window;
    if j + tdoa_window > num_point 
         segment_size = num_point - j; 
    end
    
    % add hamming window
    win_data = pcm(j:j+segment_size-1, :) .* repmat(hamming(segment_size), 1, num_channel);
    data = zeros(tdoa_window, num_channel);
    data(1:segment_size, :) = win_data;
    
    % calc gccphat tdoa
%    tdoa = gccphat(data, data(:,1));
    
    % following code do mvdr beamforming 
    % 1. calc corrvariance matrix with num_corr_var

    
    %time = tdoa / fs; 
    time = zeros(1, num_channel);
    % 2. calc w from MVDR
    w = zeros(num_channel, stft_len / 2 + 1);
    for k = 1 : stft_len /2 +1
        f = k * fs /(stft_len / 2 + 1);
        alpha = exp(-i * 2 * pi * f * time)';
        r_inv = inv(corrvar(:, :, k));
        w(:, k) = r_inv * alpha / (conj(alpha') * r_inv * alpha); % MVDR
    end
   
    %w = ones(num_channel, tdoa_window) / num_channel;
    
    % 3. sum signal
    fft_data = fft(win_data, stft_len);
    rec_signal = fft_data(1:stft_len/2+1, :) .* w';
    rec_signal = [rec_signal; flipud(rec_signal(2:end-1, :))];
    res = real(ifft(sum(rec_signal, 2)));
    res = res(1:tdoa_window);
    output(j:j+segment_size -1, :) = output(j:j+segment_size -1, :) + res.*hamming(segment_size);
    %res = real(ifft(sum(fft_data .* w', 2))) ;
   % output(j:j+segment_size-1, :) = res(1:segment_size, :) ./ hamming(segment_size);
end

sound([output; pcm(:, 1)], fs)
audiowrite('2.12.mvdr.wav', output, fs);









