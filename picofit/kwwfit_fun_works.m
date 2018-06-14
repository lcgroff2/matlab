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
  yval1=zeros([length(fitinfo.fun) length(x)]);
  for cnt=1:numfun,
    if strcmp(fitinfo.fun{cnt}.function,'exp'),
      tau=fitinfo.fun{cnt}.curval;
      tmpyval=exp(-(x-x(1))/tau); % might need transpose op here
      yval1(cnt,:)=tmpyval(:)';
    elseif strcmp(fitinfo.fun{cnt}.function,'offset'), % flat offset
      tmpyval=ones(size(x));
      yval1(cnt,:)=tmpyval(:)';
    elseif strcmp(fitinfo.fun{cnt}.function,'kww'), % stretched exponential
      tau=fitinfo.fun{cnt}.curval(1);
      beta=fitinfo.fun{cnt}.curval(2);
      tmpyval=exp(-((x-x(1))/tau).^beta);
      yval1(cnt,:)=tmpyval(:)';
    elseif strcmp(fitinfo.fun{cnt}.function,'hside'),
      tmpyval=(x>0);
      yval1(cnt,:)=tmpyval(:)';
    end  % note: later will actually add in kww function, "slope", "hside".
  end
  dx=x(2)-x(1);
  A = zeros([length(x) numfun]); % hack: extra 2 for "offset" and "slope"
  
  % for each returned function, convolve with instrument response
  % function.
  irfshift=irfshifter(T0);
  ran=1:length(x);
  for cnt=1:numfun
    %yval=filter(yval1(cnt,:),1,irfshift); % note: this is "inlining" conv
    yval=myconv(yval1(cnt,:),irfshift);
    yval=yval(ran);  % might be faster to use FFT
    yval=yval(:);
    A(:,cnt)=yval;
  end
 
 c=A\y; %least squares solution for linear coefficients (amplitudes)
 z=A*c;
 AMP=c;
 sqerror=sum(((z-y).^2));
 z=A; % return the multiple convolved functions
