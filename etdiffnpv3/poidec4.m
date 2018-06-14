
% poidec4.m

% calculate exponential decay for ensemble of particles with
% poisson-distributed quenchers. This version uses tau_kww and beta
% from the exciton diffusion model, and also differentiates between
% defect-quenchers and dyes.


%tau=3000; % lifetime without quenchers, in ps
avedye=2.992; % average number of dyes per np
avedef=2.6; % average number of defect-quenchers per np
dt=1; % 1 ps
maxnumdyes=15; % only do poisson-weighted average up to this number of dyes or defects

% will not use the variables below, instead use the direct results
% of etmulti.m
%ndyevec=[0 1 2 3 4 5]; % array of number of quenchers, must be sequential integers starting at zero
%tauvec=[3000 1500 750 400 200 100]; % tau_kww from exciton diffusion simulation
%betavec=[1 0.95 0.9 0.85 0.8 0.75]; % beta values from exciton diffusion simulation
%qevec=[0 0.5 0.75 0.8 0.85 0.9]; % quenching efficiency--should go up to 0.97 or higher
tau0=plist{1}.tau0;

% some sanity-checking, if qevec(end) is too small, then halt
if qevec(end)<0.75
    disp('quenching efficiency too small, rerun etmulti with more quenchers');
    return
end

dyeint=0:maxnumdyes; % number of dyes, for interpolation, integers starting at zero
dyedat=[0;dyevec(:)]; % start at zero dyes
qedat=[0;qevec(:)]; % need to start array at zero dyes (=zero quenching)
taudat=[tau0;taukwwvec(:)]; % need to start array with unquenched lifetime, tau0
betadat=[1;betavec(:)]; % first element should be beta=1

% for now doing simple interpolating, but should really take into account the underlying
% functional character

% hack to partially handle extrapolation--for now, just use the final value
% (hopefully the cases of high numbers of dyes/quenchers don't contribute much
% to the fluorescence anyways, and we need this for convenience)
if dyevec(end)<maxnumdyes
    dyedat(end+1)=maxnumdyes;
    qedat(end+1)=qedat(end);
    taudat(end+1)=taudat(end);
    betadat(end+1)=betadat(end);
end

% now do the interpolation
qeint=interp1(dyedat,qedat,dyeint,'cubic');
tauint=interp1(dyedat,taudat,dyeint,'cubic');
betaint=interp1(dyedat,betadat,dyeint,'cubic');


% first calculate quenching efficiency based on poisson distribution,
% using nested for-loops.
qave=0;
sumprob=0;

tpoi=0:dt:(tau0*4);
poidecay=zeros(size(tpoi));

for ndye=dyeint;
    for ndef=dyeint;
        % note: the probability of having a given bumber of dyes and a given number of defects is given by the product of the poisson distribution functions 
        prob=poisspdf(ndye,avedye)*poisspdf(ndef,avedef);
        ntot=ndye+ndef;
        if ntot>=maxnumdyes,
            ntot=maxnumdyes-1; % truncate if ntot too big (stop at qe=0.97 or so)
        end
        sumprob=sumprob+prob; % keep track of overall probability, as a cross-check
        qave=qave+prob*qeint(ntot+1); % quenching efficiency depends on the total number of quenchers
        tau_n=tauint(ntot+1);
        beta_n=betaint(ntot+1);
        % calculate weighted average decay
        poidecay=poidecay+prob*exp(-(tpoi/tau_n).^beta_n);
    end
end

fprintf(1,'for avedye = %.4g, avedef = %.4g, poisson-weighted average quenching efficiency qave = %.4g\n',avedye,avedef,qave);

if sumprob<0.9
  disp('You need a longer ndyevec, tauvec, betavec, qevec')
  return
end

% fit to kww model
fit.fnstr='A*exp(-(t/tau).^beta)'; % enter the desired function here
fit.fnxstr='t';
fit.xstr='xdat';
fit.ystr='ydat';
lbound=[0.5 100 0.3];
ubound=[1.5 4000 1.2];
% lower bound for fit parameters,
% listed in order of appearance in above formula, starting on the left 
fit.lbound=lbound; 
fit.ubound=ubound;

% do some copying/rearranging of variables from the simulation:
xdat=tpoi(1:2000);
ydat=poidecay(1:2000);

rnlfit
kwwparms=outparm;
fitkww = yfit';

% next set up for bi-exponential fit
fit.fnstr='A*(f*exp(-t/tau1)+(1-f)*exp(-t/tau2))'; % enter the desired function here
lbound=[0.5*maxy 0.1 50 1000];
ubound=[1.5*maxy 0.9 500 4000];
fit.lbound=lbound; 
fit.ubound=ubound;

rnlfit
fprintf(1,'weighted average lifetime %8.1f\n',outparm(2)*outparm(3)+(1-outparm(2))*outparm(4));
biexpparms=outparm;
fitbiexp = yfit';

fprintf(1,'plotting fit to kww and bi-exponential\n');
plot(tpoi,poidecay,xdat,fitkww,xdat,fitbiexp)
axis([0 3000 0 1])

kwwsqerr=sum((fitkww-ydat).^2);
fprintf(1,'sum of squared residuals for KWW fit %8.3e \n',kwwsqerr);
biexpsqerr=sum((fitbiexp-ydat).^2);
fprintf(1,'sum of squared residuals for bi-exponential fit %8.3e\n',biexpsqerr);