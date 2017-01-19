function [rec_signal, out, ds, w] = domvdr_flat(wav_file)
%input: 
%   wav_file: muliti-channel wav name 


%output:
%   rec_signal: single enhanced channel spectrum
%   out: single channel speech
%   ds: Delay-and-Sum Beamforming output 
%   w : MVDR weights

addpath('./utils');
[data, fs] = audioread(wav_file);
[nsamples, nchs] = size(data);


tdoa = zeros(1, nchs);

stft_len = 512;
frame_len = 400;
frame_shift = 160;
num_sil = 20; % we presume that the first num_sil is noise frames
x = enframe(data(:, 1), hamming(frame_len), frame_shift);
X = zeros(nchs, size(x, 1), stft_len); 
for ch = 1:nchs
    x = enframe(data(:, ch), hamming(frame_len), frame_shift);
    X(ch, :, :) = fft(x, stft_len, 2);
end

%% estimate the Rnn using the first num_sil
psd = outProdND(X(:, 1:num_sil, :));
Rnn = mean(psd, 3); 
output = zeros(size(data, 1), 1);
w = zeros(nchs, stft_len / 2 + 1);
i = sqrt(-1);
for k = 0 : stft_len / 2        

    f = k * fs / stft_len;
    alpha(:, k+1) = exp(-i * 2 * pi * f * tdoa');
    r_inv = inv(squeeze(Rnn(:, :, k+1)));
    w(:, k+1) = r_inv * alpha(:, k+1) / (alpha(:, k+1)' * r_inv * alpha(:, k+1)); % MVDR
    %%%%%%%%%%delay & sum%%%%%%%
    
end    

% 3. sum signal
rec_signal = zeros(stft_len/2+1, size(x, 1));
for t = 1:size(x, 1)
    Xt = squeeze(X(:, t, 1:stft_len/2+1));
    Xt = mean(Xt .* conj(alpha), 1);
    DS(t, :) = Xt;
    rec_signal(:, t) = sum(w' .* (squeeze(X(:, t, 1:stft_len/2+1))).', 2);
end    
out = OverlapAdd2(abs(rec_signal), angle(rec_signal), frame_len, frame_shift);  
ds = OverlapAdd2(abs(DS.'), angle(DS.'), frame_len, frame_shift);  



