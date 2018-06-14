function [sqerror,z,AMP]=picofit_fun(p)
  % PICOFIT_FUN.M
  %
  % function minimized by PICOMIN in the PICOFIT routine.  Performs
  % the convolution, etc., saving the user from the work.
  % The function name is carried in the global variable Fstring.
  %
  % TODO: check globals, if-then for whether to use GAUSS, OFFSET
  global IRFx IRFy FixedTaus Xdata Ydata
  x=Xdata;
  y=Ydata;
  
  T0=p(1); % T0 is the first parameter.	       
  
  taus=p(2:length(p));
  
  yval1 =feval('picoexp',x-x(1),[FixedTaus taus]);
  sz=size(yval1);
  dx=x(2)-x(1);
  numfun=sz(1); % number of vectors passed (three for triexponential)
  A = zeros([length(x) numfun+2]); % hack: extra 2 for "offset" and "slope"
  
  % for each returned function, convolve with instrument response
  % function.
  irfshift=irfshifter(T0);
  ran=1:length(x);
  for cnt=1:numfun
    yval=filter(yval1(cnt,:),1,irfshift); % note: this is "inlining" conv
    yval=yval(ran);  % might be faster to use FFT
    yval=yval(:);
    A(:,cnt)=yval;
  end
 
 A(:,numfun+1)=ones(size(yval)); % hack to deal with constant offset
 A(:,numfun+2)=(1:length(yval))'; % hack to deal with sloping offset
 c=A\y; %least squares solution for linear coefficients (amplitudes)
 z=A*c;
 AMP=c;
 sqerror=sum(((z-y).^2));
 z=A; % return the multiple functions
