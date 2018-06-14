function [sqerror,z,AMP]=kwwfit_fun(fitinfo)
  % KWWFIT_FUN.M
  %
  % function minimized by PICOMIN in the PICOFIT routine.  Performs
  % the convolution, etc., saving the user from the work.
  % The function name is carried in the global variable Fstring.
  %
  % TODO: check globals, if-then for whether to use GAUSS, OFFSET
  global IRFx IRFy FixedTaus Xdata Ydata
  x=Xdata;
  y=Ydata;
  
  T0=fitinfo.t0curval;
  
  %taus=p(2:length(p));
  % replace line below with a loop over all functions
  % also need to make a "kww" function and a simpler "exp" function
  % or just "special case" the 'exp' and the 'kww'
  numfun=length(fitinfo.fun); % number of vectors passed (three for triexponential)
  %yval1=zeros([length(fitinfo.fun) length(x)]);
  A = zeros([length(x) numfun]);
  irfshift=irfshifter(T0);

  ran=1:length(x);
  ran=ran(:);
  for cnt=1:numfun,
    if strcmp(fitinfo.fun{cnt}.function,'exp'),
      tau=fitinfo.fun{cnt}.curval;
      yval=exp(-(x-x(1))/tau);
      yval=myconv(yval,irfshift); % convolve the exponential
      yval=yval(ran);
      yval=yval(:);
      A(:,cnt)=yval;
    elseif strcmp(fitinfo.fun{cnt}.function,'offset'), % flat offset
      yval=ones(size(x));
      yval=yval(:); % just a constant offset, do not convolve
      A(:,cnt)=yval;
    elseif strcmp(fitinfo.fun{cnt}.function,'kww'), % stretched exponential
      tau=fitinfo.fun{cnt}.curval(1);
      beta=fitinfo.fun{cnt}.curval(2);
      yval=exp(-((x-x(1))/tau).^beta);
      yval=myconv(yval,irfshift); % convolve the exponential
      yval=yval(ran);
      A(:,cnt)=yval;      
    elseif strcmp(fitinfo.fun{cnt}.function,'hside'),
      cirf=cumsum(irfshift);
      cirf=cirf/cirf(end);
      cirf=cirf(:);
      A(:,cnt)=cirf;
    elseif strcmp(fitinfo.fun{cnt}.function,'irf'),
      A(:,cnt)=irfshift(:);
    end  
  end
 
  c=A\y; %least squares solution for linear coefficients (amplitudes)
  z=A*c;
  AMP=c;
  sqerror=sum(((z-y).^2));
  z=A; % return the multiple functions
