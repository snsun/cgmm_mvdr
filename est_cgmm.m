function [ lambda_v, lambda_y, Ry, Rv ] = est_cgmm( ffts )
%EST_CGMM is used to estimate the Complex GMM parameters 
%and generate the mask for nosie only and noisy t-f bins
%   ffts: M*L*(fft_len/2+1), the multi-channel fft matrix
%   lambda_v: the mask for noise only t-f bins
%   lambda_y: the mask for noisy t-f bins
%   Ry, Rv: the spacial covariance matrix of noisy and noise;;
%           M*M*F;

[M, T, F ] = size(ffts);

lambda_v = zeros(T, F);
lambda_y =zeros(T, F);
Ry = squeeze(mean(outProdND(ffts), 3));
Rv = eye(M);
Rv = reshape(Rv, [size(Rv, 1), size(Rv, 2), 1]);
Rv = repmat(Rv, [1, 1, F]);
phi_y = ones(T, F);
phi_v = ones(T, F);


for iter=1:10
    for f=1:F
        Ry_f = Ry(:, :, f);
        Rv_f = Rv(:, :, f);
        y_tf = ffts(:, :, f);
        y_y_tf = outProdND(y_tf);
        sum_y = zeros(M);
        sum_v = zeros(M);
        e= eye(M)*0.00001;
        for t = 1:T
            phi_y(t, f) = 1/M*trace(y_y_tf(:, :, t)*inv(Ry_f+e));
            phi_v(t, f) = 1/M*trace(y_y_tf(:, :, t)*inv(Rv_f+e));    
            kernel_y = y_tf(:, t)' * inv(phi_y(t, f)*Ry_f+e) * y_tf(:, t);
            kernel_v = y_tf(:, t)' * inv(phi_v(t, f)*Rv_f+e) * y_tf(:, t);
            p_y(t, f) = exp(-kernel_y)/(pi*(det(phi_y(t, f)*Ry_f)));
            p_v(t, f) = exp(-kernel_v)/(pi*det(phi_v(t, f)*Rv_f));
            lambda_y(t, f) = p_y(t, f) / (p_y(t, f)+p_v(t, f));
            lambda_v(t, f) = p_v(t, f) / (p_y(t, f)+p_v(t, f));
            sum_y = sum_y + lambda_y(t, f)/phi_y(t, f)*y_y_tf(:, :, t);
            sum_v = sum_v + lambda_v(t, f)/phi_v(t, f)*y_y_tf(:, :, t);
        end
        Ry(:, :, f) = 1/sum(lambda_y(:, f)) * sum_y;
        Rv(:, :, f) = 1/sum(lambda_v(:, f)) * sum_v;
        
        
        
    end
    Q = sum(sum(lambda_y .* log(p_y+0.001) + lambda_v .* log(p_v+0.001)))
    figure(1)
    imagesc(real([flipud(lambda_y');flipud(lambda_v')]));
end  


end

