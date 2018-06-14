function y=jmgauss(xtmp,fwhm)
%
%jmgauss(x,fwhm): a regular gaussian.
%
% x is the x vector, sig is "sigma", FHWM = 2.355*sigma
%
sig=fwhm/2.355;
y=exp(.5/(sig*sig)*(-(xtmp.*xtmp)))/sig*sqrt(2*pi)/2/pi;
