
% script to perform fitting on energy transfer model simulation
% results. Uses RNLFIT, and data from simulation.


% do some copying/rearranging of variables from the simulation:
xdat=tt;
ydat=fdecay;
maxy=max(ydat);
xdat=xdat(2:end); % first point is bad due to the way histogram is calculated
ydat=ydat(2:end);

% first set up for single-exponential fit
fit.fnstr='A*exp(-t/tau)'; % enter the desired function here
fit.fnxstr='t';
fit.xstr='xdat';
fit.ystr='ydat';
lbound=[0.5*maxy 50];
ubound=[1.5*maxy 4000];
% lower bound for fit parameters,
% listed in order of appearance in above formula, starting on the left 
fit.lbound=lbound; 
fit.ubound=ubound;



if exist('parm','var')
  clear npparm
  npparm=parm; % this script uses "parm" in a different way -- fit parameters,
             % so we'll move them them into variable "npparm"
  clear parm
end
rnlfit
expparms=outparm;
if exist('plist','var'),
  tauvec(pidx)=outparm(2);
end


%return
if exist('plist','var')==0
  disp('press any key for next fit')
  pause
end

% next set up for bi-exponential fit
fit.fnstr='A*(f*exp(-t/tau1)+(1-f)*exp(-t/tau2))'; % enter the desired function here
lbound=[0.5*maxy 0.1 50 1000];
ubound=[1.5*maxy 0.9 500 4000];
fit.lbound=lbound; 
fit.ubound=ubound;

rnlfit
fprintf(1,'weighted average lifetime %8.1f\n',outparm(2)*outparm(3)+(1-outparm(2))*outparm(4));
biexpparms=outparm;
tauavg=outparm(2)*outparm(3)+(1-outparm(2))*outparm(4);
if exist('plist','var'),
  tau1vec(pidx)=outparm(3);
  tau2vec(pidx)=outparm(4);
  fvec(pidx)=outparm(2);
  tauavgvec(pidx)=tauavg;
end

if exist('plist','var')==0
  disp('press any key for next fit')
  pause
end

fit.fnstr='A*exp(-(t/tau).^beta)';
lbound=[0.5*maxy 50 0.3];
ubound=[1.5*maxy 900 1.1];
fit.lbound=lbound; 
fit.ubound=ubound;
rnlfit;
strexpparms=outparm;

if exist('plist','var'),
  taukwwvec(pidx)=outparm(2);
  betavec(pidx)=outparm(3);
  qevec(pidx)=mean(qyet);
end
