function yout=nac(y)

% function for "properly normalized" autocorrelation, for
% obtaining number of molecules in excitation volume.

% returns "yout".  In principle, yout(1) is the reciprocal of
% the mean number of fluorophores.  However, Poisson and other 
% noise also affect yout(1), so it is best to estimate G(0) by
% fitting the data points *after* yout(1).

yout=acft(y)/mean(y)^2/length(y)-1;

