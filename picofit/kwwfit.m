function [bestfitinfo,yfit,AMP]=kwwfit(x,y,fitinfo)
% KWWFIT Nonlinear curve fitting and plotting routine with convolutions
%        Modified for stretched exponentials

% Usage: [pout,yfit,amp]=kwwfit(x,y,[T0min T0max],fitinfo);
%
%       KWWFIT(X,Y,T0,fitinfo) fits functions to a data
%       set with convolutions.
%       It uses random number generation
%       for nonlinear components and 
%       linear least squares for linear components (amplitudes).
%
%       It is assumed that x starts at zero for purposes of
%       evaluating the fit function.
%
%       The global variables IRFx,IRFy should be the instrument
%       response function.
%
%
%	function also returns YFIT, the fit curve, if required.
%
%
%       Assume the instrument function is in the global variables
%       IRFx, IRFy
%
%	Example:
%	[pfm,yfit,AMP]=picofit(DataX,DataY,[-20 20],[],[1 8 50]);
%
%	Jason McNeill, Oct 2001 (mcfit). August 2006 (picofit)
%       May 2008 (kwwfit)

%	TODO: -remove kludges and kruft
%             -check weights code
%	      -add better search method?
%             -use an object tree interface? (fitinfo)
%             -be careful about the function starting at x=0 (checked--OK)

% fitinfo is a struct with an embedded list like:
% fitinfo.t0.curval
% fitinfo.t0.lbound
% fitinfo.t0.ubound
% fitinfo.fun{1}.function = 'kww'
% fitinfo.fun{1}.lbound = [0.5 10]
% fitinfo.fun{1}.ubound = [1.1 50]
% fitinfo.fun{1}.curval = [ 0.7 30] -- current guess

% for "fixed" parameters, we can just set the lower and upper bounds
% equal (need to make sure the code handles that properly).

global IRFx IRFy Xdata Ydata Weight

IRFy=IRFy(:)';

Xdata=x(:);
Ydata=y(:);

%Weight=(abs(Ydata)+3).^(-1); % *** Seems wrong ***
%Weight=Weight/sum(Weight)*length(Weight);
Weight=ones(size(Xdata));

%taus variable no longer needed
%T0=fitinfo.t0.curval;

disp('starting...');
if nargin~=3, 
   error('usage: kwwfit(x,y,fitinfo)'); 
end
x=x(:);
y=y(:);

if length(IRFx) < 5,
  disp('need to define instrument function first and make it global')
  stop
end

% set the initial guess to be in the middle of the range:
justplot=0; % if range is zero, then just plot and return
%for cnt=1:length(fitinfo.fun)
%  fitinfo.fun{cnt}.curval = (fitinfo.fun{cnt}.lbound+fitinfo.fun{cnt}.ubound)/2;
%  if fitinfo.fun{cnt}.curval ~= fitinfo.fun{cnt}.lbound,
%    justplot=0;
%  end
%end

% now to plot the initial guess
[sqerror,yfit,AMP]=kwwfit_fun(fitinfo);

bestfitinfo=fitinfo;

figure(1)
plot(x,y,'o',x,yfit*AMP,'g',x,yfit*AMP-y(:))
%chisq=sum(((y(:)-yfit*AMP).^2).*Weight)/length(y);
%title(['mean square error =' num2str(chisq)])
ts=sprintf('RMS resid %5.2e',sqrt(sum(((y(:)-yfit*AMP).^2).*Weight))/length(y));
title(ts)
drawnow

% if there are no nonlinear parameters to adjust, return.
if justplot==1,
  %disp('just plotting, not fitting, for this call')
  return
end

%p_in=[T0 taus]

newfitinfo = kwwmin(bestfitinfo);

[sqerror,yfit,AMP]=kwwfit_fun(newfitinfo);

bestfitinfo=newfitinfo; 
   
% plot intermediate results

plot(x,y,'o',x,yfit*AMP,'g',x,yfit*AMP-y(:))
ts=sprintf('RMS resid %5.2e',sqrt(sum(((y(:)-yfit*AMP).^2).*Weight))/length(y));
title(ts)


% return



