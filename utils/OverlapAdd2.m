function ReconstructedSignal=OverlapAdd2(XNEW,yphase,windowLen,ShiftLen);

%Y=OverlapAdd(X,A,W,S);
%Y is the signal reconstructed signal from its spectrogram. X is a matrix
%with each column being the fft of a segment of signal. A is the phase
%angle of the spectrum which should have the same dimension as X. if it is
%not given the phase angle of X is used which in the case of real values is
%zero (assuming that its the magnitude). W is the window length of time
%domain segments if not given the length is assumed to be twice as long as
%fft window length. S is the shift length of the segmentation process ( for
%example in the case of non overlapping signals it is equal to W and in the
%case of %50 overlap is equal to W/2. if not givven W/2 is used. Y is the
%reconstructed time domain signal.
%Sep-04
%Esfandiar Zavarehei

if nargin<2
    yphase=angle(XNEW);
end
if nargin<3
    windowLen=size(XNEW,1)*2;
end
if nargin<4
    ShiftLen=windowLen/2;
end
if fix(ShiftLen)~=ShiftLen
    ShiftLen=fix(ShiftLen);
    disp('The shift length have to be an integer as it is the number of samples.')
    disp(['shift length is fixed to ' num2str(ShiftLen)])
end

[FreqRes FrameNum]=size(XNEW);

Spec=XNEW.*exp(j*yphase);

if mod(windowLen,2) %if FreqResol is odd
    Spec=[Spec;flipud(conj(Spec(2:end,:)))];
else
    Spec=[Spec;flipud(conj(Spec(2:end-1,:)))];
end
sig=zeros((FrameNum-1)*ShiftLen+windowLen,1);
weight=sig;
for i=1:FrameNum
    start=(i-1)*ShiftLen+1;
    spec=Spec(:,i);
    ifft_tmp = real(ifft(spec, 512));
    win = hamming(400);
%      win2=win.^2;
%  w0=mean(sum(reshape(win2,200,2),2)).^0.5;
% 
%  win=win/w0;
    sig(start:start+windowLen-1)=sig(start:start+windowLen-1)+ifft_tmp(1:400).*win;
end
ReconstructedSignal=sig;

end

