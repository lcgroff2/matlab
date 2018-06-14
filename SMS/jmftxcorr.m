function [c] = acft(x,y)

% autocorrelation using Fourier transform,
% since the autocorr. is the FT of the power spectrum.

%x=x(:)';
%n = length(x);
%y = x(n:-1:1);
%y = [y y];

%c = conv(x, y);
%c = c(n:1.5*n-1);

x=x(:)';
lnx=length(x);
%padlen=floor(lnx/2);  % better agreement with Dehong's routine
                       % if we don't pad

% no padding and FFT implies periodic boundary conditions
xft=fft(x);
yft=fft(y);
c=real(ifft(xft.*conj(yft)));
