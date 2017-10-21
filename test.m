%% Test scripts 
%% Author Sining Sun (NWPU)
% snsun@nwpu-aslp.org
clc
clear

%% Load the test multi-channel test data
M = 6; %channels number
for i = 1:M
    wav(:, i) = audioread(['test_wav/test2/F01_050C0103_STR.CH' int2str(i) '.wav']);
end

%% enframe and do fft
frame_length = 512;
frame_shift = 256;
fft_len = 512;
[frames, ffts] = multi_fft(wav, frame_length, frame_shift, fft_len);

%% Estimate the TF-SPP and spacial covariance matrix for noisy speech and noise 
[lambda_v, lambda_y, Ry, Rv] = est_cgmm(ffts);

Rx = Ry -0.8*Rv;         %trade off. Rx may be not positive definite

[M, T, F]  = size(ffts); %fft bins number
d = zeros(M, F);         %steering vectors
w = d;                   %mvdr beamforming weight 
output = zeros(T, F);    %beamforming outputs
e = 0.0001*eye(M);       %avoid the matrix singular

%% Get steering vectors d using eigvalue composition 
for f= 1:F
    [V, D] = eig(squeeze(Rx(:, :, f)));
    d(:, f) = V(:, 1);
end
%% Do MVDR beamforming
output = mvdr(ffts, Rv, d);

%% Reconstruct time domain signal using overlap and add;
output = [output, fliplr(conj(output(:, 2:end-1)))];
rec_frames = real(ifft(output, fft_len, 2));
sig = overlapadd(rec_frames, hamming(frame_length, 'periodic'), frame_shift);

audiowrite('output.wav', sig, 16000);