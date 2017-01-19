function  B = beampattern( w)
%BEAMPATTERN Summary of this function goes here
%   Detailed explanation goes here

[x, y, z] = sphere(50);
x = x;
y = y - 0.095;
z = z ;
% xmic=[-0.10 0 0.10]'; % left to right axis
% ymic=[-0.095 -0.095 -0.095]'; % bottom to top axis
% zmic=[ 0 0 0]'; % back to front axis

xmic=[-0.10  0.10 -0.10 0 0.10]'; % left to right axis
ymic=[0.095  0.095 -0.095 -0.095 -0.095]'; % bottom to top axis
zmic=[0  0 0 0 0]'; % back to front axis

dd = ones(size(h));
fs = 16000;
stft_len = 512;
for m = 1:size(x, 1);
    for n = 1:size(x, 2)  
        
        for i = 1:size(zmic, 1)
           real_tdoa(i) = norm([x(m, n), y(m, n), z(m, n)] - [xmic(i), ymic(i), zmic(i)]);
        end
        real_tdoa(:) = (real_tdoa(:) - real_tdoa(1))/334;
        
%         ang = [x(m, n), y(m, n), z(m, n)]*[xmic(end), 0 ,zmic(end)]'...
%             /(norm([x(m, n), y(m, n), z(m, n)]) * norm([xmic(end), 0 ,zmic(end)]));
%         real_tdoa = [0, 1, 2]*0.1*ang/343;
% %           
        k = 1:stft_len/2+1;
        f = (k-1)' * fs / stft_len;
        d = repmat(f, [1, size(zmic, 1)]) .* repmat(real_tdoa, [stft_len/2+1, 1]); %257*5
        d = exp(-(-1)^(0.5) * 2 * pi * d'); %5*257
        %B(m, n, :) = sum(conj(d) .* h, 1);
        B(m, n, :) = sum(conj(d) .* h, 1);
        [TH(m, n), PHI(m, n), R(m, n)] = cart2sph(x(m, n), y(m, n), z(m, n));
        %[Bx(m,n), By(m, n), Bz(m, n) ]= sph2cart(TH(m, n), PHI(m, n), abs(B(m, n, K)));
        Bz(m,n) = abs(B(m, n, K))*sin(PHI(m,n));
        Bx(m, n) = abs(B(m, n, K))*cos(PHI(m, n))*cos(TH(m, n));
        By(m, n) = abs(B(m, n, K))*cos(PHI(m, n))*sin(TH(m, n));
        
        %[Bx(m,n), By(m, n), Bz(m, n) ]= sph2cart(TH(m, n), PHI(m, n), R(m, n));
    end
end
mesh(Bx, By, Bz);


x = [0, p(1)-3]; y = [0, p(2)-2.5]; z = [0, p(3)-1.5];
line(1.5*x, 1.5*y, 1.5*z);

title(num2str(K*fs/512));

end

