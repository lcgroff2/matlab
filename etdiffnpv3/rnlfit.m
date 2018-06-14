
% RNLFIT.M -- "random nonlinear fit".  Generate random numbers
% to use as initial values for a nonlinear fit.

% also has bits to help try to detect local minima
% and other places where the fit routine might get hung up,
% by generating a long list of random initial parameters,
% feeding them to the nonlinear fit routine, and keeping
% track of where the fits "end up" in each case, along with
% the square error.  Then can generate a histogram of square
% error values and parameters, to see if there are multiple
% "peaks" in the distribution.

% this script depends on helper script: fit_parsefun

% assumes x (or r) data in variable "xdat", y data in variable
% "ydat", and function and other parameters defined in a struct like:
%fit.fnstr='A*exp(-t/tau)'; % enter the desired function here
%fit.fnxstr='t';
%fit.xstr='xdat';
%fit.ystr='ydat';
%fit.lbound=[200 300];  
%fit.ubound=[1200 4000];

% This script calculates the RMS residual (square error) per data point,
% to make it easier to compare to what is expected, based on poisson
% statistics.


fit.iniparm=(fit.lbound+fit.ubound)/2; % initial values for parameters
nguess=200; % number of initial guesses.

% allocate arrays for storing fit results, later histogramming
parmsmat=zeros([nguess length(fit.lbound)]);
%fsparmsmat=parmsmat;
sqerrvec=zeros([nguess 1]); % also save sqerr
%fssqerrvec=sqerrvec;
errcnt=0; % keep track of how many failed guesses
badparms=[];
%fsbadparms=[];

if ( ~exist('xdat') ) || ( ~exist('ydat') )
  disp('you need to enter xdat and ydat first');
  return
end

fprintf(1,'function is: %s\ntrying to parse...',fit.fnstr);

% "Parse" function, replacing parameters with parm(1), parm(2), ...
fit=fit_parsefun(fit);


% plot data
axis normal % automatic axis adjustment
eval(['plot(' fit.xstr ',' fit.ystr ',''o'')']);
hold on

% basic check of parameter inputs
if length(fit.parmlist) ~= length(fit.iniparm)
  fprintf(1,'error: the number of parameters does not match the number of initial values\n');
  fit.parmlist
  fit.lbound
  return;
end

% repeat back the parameters and their initial values
fprintf(1,'function parsed\nparameters and initial values:\n');
for idx=1:length(fit.parmlist)
  fprintf(1,'%s = %f \n', fit.parmlist{idx},fit.iniparm(idx));
end


% next, set up anonymous function for fitting, and try to 
% evaluate the function and plot it
inparm=fit.iniparm;
xvar=eval(fit.xstr);
xvar=xvar(:); % nlinfit etc. require column vectors, maybe?
yvar=eval(fit.ystr);
yvar=yvar(:);

myfunstr=['@(parm,' fit.xstr ')(' fit.parsedfnstr ')'];
myfun=eval(myfunstr);
yfit=myfun(fit.iniparm,xvar);
if length(yfit) ~= length(xvar),
  disp('problem evaluating function, edit function and try again');
  return
end

% set up anomymous function for fsolve
%fsfunstr=['@(parm)(' fit.parsedfnstr '-' fit.ystr ')'];
%fsfun=eval(fsfunstr);

% set up anonymous function for fminsearch
%fmfunstr=['@(parm)(sum((' fit.parsedfnstr '-' fit.ystr ').^2))'];
%fmfun=eval(fmfunstr);

if sum(isnan(yfit)+isinf(yfit))>0, % if bad pre-fitting results
  disp('Evaluating function with the initial parameters given');
  disp('yields infinity...edit parameters and try again');
  return
end


% next ask if guesses look OK, if so press any key else ctrl-c
plot(xvar,yfit)
hold off
%fprintf(1,'Check plot to make sure guesses look OK.  If OK press any key,\n');
%fprintf(1,'else press Ctrl-c and edit initial guess parameters\n');
%pause

isoct=exist('OCTAVE_VERSION','builtin');

% start loop of random initial guesses
for ridx=1:nguess
  
  % generate random initial guess parameters
  inparm=fit.lbound+rand(size(fit.lbound)).*(fit.ubound-fit.lbound);
  
  % now run nlinfit in matlab, or nonlin_curvefit in Octave
  if isoct==0,
    try
      [outparm,resid,J,Sigma]=nlinfit(xvar,yvar,myfun,inparm);
    catch
      resid=Inf;
      outparm=inparm;
      errcnt=errcnt+1;
    end
  else
    try
      optset=optimset('TolFun',1e-13);
      [outparm,fity,cvg,numits]=nonlin_curvefit(myfun,inparm(:),xvar,yvar,optset);
      resid=yvar(:)-fity(:);
      %fprintf(1,'number of iterations: %i\n',numits.niter);
    catch
      %disp('error in nonlincurvefit for guess params:')
      %inparm
      resid=Inf;
      errcnt=errcnt+1;
      outparm=inparm;
    end
  end
  sqerr=sum(resid.*resid);
  if (imag(sqerr)^2>1e-16*real(sqerr)^2) || (imag(sqerr)^2>1e-30)
    %disp('nlinfit yielded a complex square error for guess params:')
    %inparm
    resid=Inf;
    sqerr=Inf;
    outparm=inparm;
    errcnt=errcnt+1;
  end
  


  % now try with fsolve - which doesn't handle xdat, ydat as 
  % nlinfit does, so we need a different type of anonymous function

  %try
    %warn off % matlab displays some useless warnings
  %  [fsoutparm,fval,info,output]=fsolve(fsfun,inparm);
    %warn on
  %  parm=fsoutparm;
  %  fsfity=eval([fit.parsedfnstr ';']);
  %  fsresid=yvar(:)-fsfity(:);
  %catch
  %  fsresid=Inf;
    %disp('fsolve terminated in an error');
  %  fsoutparm=inparm;
  %  errcnt=errcnt+1;
  %end
  %fssqerr=sum(fsresid.*fsresid);
  %if (imag(fssqerr)^2>1e-16*real(fssqerr)^2) || (imag(fssqerr)^2>1e-30)
    %disp('fsolve yielded a complex square error for guess params:')
    %inparm
  %  fsresid=Inf;
  %  fssqerr=Inf;
  %  fsoutparm=inparm;
  %  errcnt=errcnt+1;
  %end

  %if (abs(fssqerr-sqerr)>0.01*sqerr) && (isinf(fssqerr)+isinf(sqerr)<1),
    %disp('WARNING: fsolve and nlinfit/nonlin_curvefit yielded different results.  Check parameters');
  %end

  sqerrvec(ridx)=sqerr;
  %fssqerrvec(ridx)=fssqerr;
  parmsmat(ridx,:)=outparm;
  %fsparmsmat(ridx,:)=fsoutparm;
end
  
% determine yfit
[minval,minidx]=min(sqerrvec);
outparm=parmsmat(minidx,:);
sqerr=minval;

%[fsminval,fsminidx]=min(fssqerrvec);
%if fsminval<minval
%  outparm=fsparmsmat(fsminidx,:);
%  sqerr=fsminval;
%end

yfit=myfun(outparm,xvar);

% first clean up parmsmat and fsparmsmat, removing bad points
badidx=find(isinf(sqerrvec));
nbad=length(badidx);
if length(badidx)>0
  for idx=1:length(badidx)
    badparms=[badparms parmsmat(badidx(idx),:)]; % make a record
    parmsmat(badidx(idx),:)=outparm; % replace with good values
  end
end
parmsmat=real(parmsmat); % nlinfit sometimes gives an imaginary part
%badidx=find(isinf(fssqerrvec));
%fsnbad=length(badidx);
%if length(badidx)>0
%  for idx=1:length(badidx)
%    fsbadparms=[fsbadparms fsparmsmat(badidx(idx),:)];
%    fsparmsmat(badidx(idx),:)=outparm; % replace with good values
%  end
%end
    
%disp('currently plotting initial guess of middle of each value range...')
disp('');
disp('the stdev reported below is not the error bars, but rather a');
disp('measure of how often the fit routine converges to the same point.');
disp('ideally, the stdev below should be zero. If std is more than');
disp('1e5*param, then there are likely multiple minimima.');
for idx=1:length(fit.lbound),
  %disp('press any key for (next) histogram')
  %pause
  %[NN,XX]=hist(parmsmat(:,idx));
  %plot(XX,NN)
  fprintf(1,'For parameter "%s", mean fit value is %8.5g, stdev in fit values is %7.4g\n', fit.parmlist{idx},mean(parmsmat(:,idx)),std(parmsmat(:,idx)));

end

% plot fit, keeping same axes
hold off
axis normal % hold axes
plot(xvar,yvar,'o')
hold on
ax=axis;
plot(xvar,yfit)
axis(ax);
hold off

% Now display final parameter values 
fprintf(1,'Best Fit results...\n')
for idx=1:length(outparm)
  fprintf(1,'%s = %8.5g\n',fit.parmlist{idx},outparm(idx));
end
fprintf(1,'RMS residual is %7.4g per data point, in y axis units (e.g., counts).\n',sqrt(sqerr/length(yfit)));
fprintf(1,'The underlying poisson noise (theoretical residual) is estimated at %7.4g counts per data point\n',mean(sqrt(yfit)));
fprintf(1,'Relative RMS residual ( rms/<data> ) is %4.1f percent\n',100*sqrt(sqerr/length(yfit))/mean(ydat));
fprintf(1,'nlinfit/nonlinfit: %5.2f pct of random initial guesses failed to converge\n',nbad/nguess*100);
%fprintf(1,'fsolve: %5.1f pct of random initial guesses failed to converge\n',fsnbad/nguess*100);

% reshape bad parms vector into a matrix, for possible later plotting
% to determine the "badlands" of parameter space
bpmat=reshape(badparms,length(fit.lbound),[]);

% plot(bpmat(3,:),bpmat(4,:))