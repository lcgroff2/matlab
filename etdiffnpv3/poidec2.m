% calculate exponential decay for ensemble of particles with
% poisson-distributed quenchers. This version uses tau_kww and beta
% from the exciton diffusion model




tau=800; % lifetime without quenchers
nbar=2.992; % average number of quenchers per np

ndyevec=[0 1 2 3 4 5]; % array of number of quenchers, must be integers
tauvec=[800 571 324 230 153 111]; % tau_kww from exciton diffusion simulation
betavec=[1 .7517 .6893 .6820 .6539 .6185]; % beta values from exciton diffusion simulation
qevec=[0 .7090 .8826 .9269 .9594 .9746]; % quenching efficiency

% first calculate quenching efficiency based on poisson distribution
poidist=poisspdf(ndyevec,nbar);
qave=sum(poidist.*qevec);
fprintf(1,'for nbar = %f, poisson-weighted average quenching efficiency qave = %f\n',nbar,qave);

if sum(poidist)<0.9
  disp('You need a longer ndyevec, tauvec, betavec, qevec')
  return
end


tpoi=0:2:(tau*4);
poidecay=zeros(size(tpoi));

for idx=1:length(ndyevec), 
  n=ndyevec(idx);
  p_n=poisspdf(n,nbar);
  %f_n=(1-qevec(idx)); % fluorescence, relative to n=0, for n quenchers
  %q_ave=q_ave+(1-f_n)*p_n; % calculate overall ave quenching efficiency
  %tau_n=tau*f_n; % I believe this is correct, for lifetime in presence of quencher
  %tauvec(idx)=tau_n;
  tau_n=tauvec(idx);
  beta_n=betavec(idx);
  poidecay=poidecay+p_n*exp(-(tpoi/tau_n).^beta_n);
end

% fit to kww model
fit.fnstr='A*exp(-(t/tau).^beta)'; % enter the desired function here
fit.fnxstr='t';
fit.xstr='xdat';
fit.ystr='ydat';
lbound=[0.5 500 0.3];
ubound=[1.5 4000 1.2];
% lower bound for fit parameters,
% listed in order of appearance in above formula, starting on the left 
fit.lbound=lbound; 
fit.ubound=ubound;

% do some copying/rearranging of variables from the simulation:
xdat=tpoi(1:1600);
ydat=poidecay(1:1600);

rnlfit
kwwparms=outparm;

fprintf(1,'plotting fit to kww\n');
plot(tpoi,poidecay,xdat,yfit)
axis([0 3000 0 1])




