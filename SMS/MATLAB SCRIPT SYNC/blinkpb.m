
% blinktst - new script for simulating blinking. I decided to start a new script rather 
% than edit "blinkdemo". This script is for testing several ideas and/or answering questions:
% -whether it makes a difference if FT is a power of two or not.
% -whether quenching efficiency of "polaron" affects C(tau)/autocorr/acft/acvec unnormalized correlation, or xcorr or nac G(tau) normalized correlation
% -effect of offset, readout noise, poisson noise, and ADC gain on C(tau) or G(tau)
% -whether C(1) is the sum of all data squared, for CMOS data (we know this is true for APD data)
% -what is the minimum signal level and number of frames to get meaningful autocorrelation results
% -try to use 2d autocorrelation or other approach to separate the two rate constants
% -check the effects of weak blinking (set QE=0.05 or similar) on autocorrelation, etc.
% -check to see if quenching efficiency can be extracted from the autocorrelation using approach in nano lett paper.
% -test photon counting histogram approach for finding equilibrium constant (and separate out rates)
% -test whether Barbara approach to finding equilibrium constant will work for CMOS data
% -test histogram of "on" times and "off" times
% -add photobleaching to script

% Note1: this version does some "oversampling" in time, to approximate the effect of a molecule switching from one state to another midway through a time step. This also allows us to examine the effect of having one or both rates "fast" relative to framerate. The oversampling is moderate, (10x), but that should be enough to give an approximately correct flavor to the data.

% Note2: simulated cmos data is in the variable "simcmos". It is simulating "one pixel" of scmos camera, which is of course different than what you get from the tracking code, which fits a range of pixels to a gaussian.

% Note3: I decided to put the plotting and autocorrelation in "blinktstnotes.txt"

% TODO: add code for those tests

% future: need to do similar tests of simulated MSD, or even simulated CMOS spots

% experimental parameters
dt = 1e-3;

tmax=120; % number of seconds for experiment

%tau_on = 0.010; % lifetime of "on" state

tau_off = 0.005; % lifetime of "off" state

emrate=5000000; % number of photons per second, when in the "on" state

%QE = 0.1; % quenching efficiency ( around 1.0 for single dye )


% camera settings and properties
offset=0;
readnoise=3;
sCMOS_quanteff = 0.57; % This is the quantum efficiency of the sCMOS camera
collect_eff = 0.03; % This is the collecting efficiency of the microscope
gain=1.6; % gain of 0.3 means roughly 3 photons per ADC count

% simulation parameters (affect how things are calculated, but do not correspond to physical picture)
oversamp=10; % basically, uses a much smaller dt when calculating the trajectory, then samples back down to dt

% derived parameters, based on the above. Do not edit.
t=dt:dt:tmax; % time array

kon = 1 / tau_off; % kon = rate const for switching on

koff = 1/ tau_on ; % koff = rate const for switching off

Keq_on = kon / koff; % equilibrium constant, [on] / [off] or n_on / n_off


ktot = kon+koff; % total rate constant

nblink_est = tmax * ktot / 4; % estimated number of blink (on+off) events

photobleach = 0;
paraout = [];
for loop=1:100
    rng('shuffle')
    % use exponentially-distributed random numbers to generate on and off times
    % should make a few extra 
    ontimes = exprnd(tau_on,ceil(nblink_est*1.2),1);
    offtimes = exprnd(tau_off, ceil(nblink_est*1.2),1);

    % just for debugging
    %fprintf('total amount of time in on and off states: %.4g\n',sum(ontimes+offtimes));

    % now generate molecular state, using time oversampling 
    curt=0.0; % current time
    curidx=1; % current idx for time vector
    bidx=1; % which number of "blink" event this is
    lent=length(t);
    ts=(dt/oversamp):(dt/oversamp):tmax;
    osstate=zeros(size(ts));


    while curt<tmax
      % assume it starts in "on" state"
      nsteps=round(ontimes(bidx)/dt*oversamp);
      curt=curt+ontimes(bidx);
      if nsteps>0 % it is possible for state to be on a very short time, corresponding to "zero" steps
        if curidx+nsteps>length(ts)
          osstate(curidx:end)=1;
        else
          osstate(curidx:(curidx+nsteps))=1;
        end
      end
      curidx=curidx+nsteps;
      % next it switches "off"
      nsteps=round(offtimes(bidx)/dt*oversamp);
      curt=curt+offtimes(bidx);
      curidx=curidx+nsteps;
      bidx=bidx+1;
      if bidx>length(ontimes),
        break;
      end
    end

    % resample molstate down so time steps are "dt"
    tmp=0;
    for idx=1:oversamp,
      tmp=tmp+osstate(idx:oversamp:end); % need to check this
    end
    molstate=tmp/10;


    % Now generate simulated data

    photonsperframe=emrate*dt;
    % use approximation that for N photons average, replace with a gaussian random number with mu=N and sigma=N
    % also take into account qe
    onpart=mypoissrnd(QE*photonsperframe*molstate); % simulated additional photons when in "on" state
    offpart=mypoissrnd((1-QE)*photonsperframe*ones(size(t))); % simulated photons when in "off" state and part of photons when in "on" state
    simphotons=round(onpart+offpart); % round to an integer number of photons

    simcmos = offset+floor(gain*simphotons*sCMOS_quanteff*collect_eff)+round(readnoise*(randn(size(t))));

    if photobleach == 1;
        % photobleaching
        pbl = rand*2*exp(-t/randi([5,200]))+1;
        simcmos = simcmos.*pbl;
    end

    %figure(1)
    %subplot(3,2,[1 2])
    %plot(t,simcmos)
    
    ran = 2:200;
    
    % 1 nac
    %acc = acft(simcmos(1:20000)-mean(simcmos(1:20000)));
    %acc = xcorr((simcmos(1:20000)-mean(simcmos(1:20000))),'coeff'); acc=acc(20000+1:end);
    acc = nac(simcmos);
    %acc = xcorr((simcmos(1:20000)-offset),'coeff');acc=acc(20000+1:end);
    acc = acc(ran);
    %acc = acc/max(acc);
    x = 1:length(acc);
    if photobleach == 1
        [f2,f2v] = fit(t(x)',acc','A*exp(-x/(B))+C*exp(-x/(D))+y0','startpoint',[0.001,0.005,0.001,1,0.001],'upper',[1,0.1,1e6,1e6,1e6],'lower',[1e-5,1e-5,-1e6,-1e6,-1e6],'tolfun',1e-9);
    elseif photobleach == 0
        [f2,f2v] = fit(t(x)',acc','A*exp(-x/(B))+y0','startpoint',[1,0.005,1],'tolfun',1e-9);
    end
    %subplot(3,2,[3 4])
    %plot(f2,t(x),acc)
    warning('off')
    fitIac = struct(f2); fitIacv = struct(f2v);
    warning('on')
    Anac = fitIac.coeffValues{1,1};
    Tnac = fitIac.coeffValues{1,2};
    if photobleach ==1
        Anac2 =fitIac.coeffValues{1,3};
    end
    Y0nac = fitIac.coeffValues{1,end};
    R2nac = fitIacv.rsquare;
    %fprintf('Ttotal(theory)=%1.2fms\n',Ttotaltheory)
    Ttotalsimnac = Tnac*1000;
%     P=polyfit(t(ran),log(nacc(ran)),1); % semilog-plot analysis, only valid if nacc goes to zero at long times, and also only if nacc(ran) is above zero everywhere
%     est_tau=-1/P(1);
    
    % 2 acft
    acc = acft(simcmos(1:20000)-mean(simcmos(1:20000)));
    acc = acc(ran);
    acc = acc/max(acc);
    x = 1:length(acc);
    if photobleach == 1
        [f2,f2v] = fit(t(x)',acc','A*exp(-x/(B))+C*exp(-x/(D))+y0','startpoint',[0.1,0.005,100,100,0],'upper',[1,0.1,1e6,1e6,1e6],'lower',[1e-5,1e-5,-1e6,-1e6,-1e6],'tolfun',1e-9);
    elseif photobleach == 0
        [f2,f2v] = fit(t(x)',acc','A*exp(-x/(B))+y0','startpoint',[1,0.005,0],'tolfun',1e-9);
    end
    %subplot(3,2,[3 4])
    %plot(f2,t(x),acc)
    warning('off')
    fitIac = struct(f2); fitIacv = struct(f2v);
    warning('on')
    Aacft = fitIac.coeffValues{1,1};
    Tacft = fitIac.coeffValues{1,2};
    if photobleach ==1
        Aacft2 = fitIac.coeffValues{1,3};
    end
    Y0acft = fitIac.coeffValues{1,end};
    R2acft = fitIacv.rsquare;
    %fprintf('Ttotal(theory)=%1.2fms\n',Ttotaltheory)
    Ttotalsimacft = Tacft*1000;
    
%     % 3 xcorr
%     acc = xcorr((simcmos(1:20000)-mean(simcmos(1:20000))),'coeff'); acc=acc(20000+1:end);
%     acc = acc(ran);
%     acc = acc/max(acc);
%     x = 1:length(acc);
%     if photobleach == 1
%         [f2,f2v] = fit(t(x)',acc','A*exp(-x/(B))+C*exp(-x/(D))+y0','startpoint',[0.1,0.005,100,1000,-100],'upper',[1,0.1,1e6,1e6,1e6],'lower',[1e-5,1e-5,-1e6,-1e6,-1e6],'tolfun',1e-9);
%     elseif photobleach == 0
%         [f2,f2v] = fit(t(x)',acc','A*exp(-x/(B))+y0','startpoint',[1,0.005,0],'tolfun',1e-9);
%     end
%     %subplot(3,2,[3 4])
%     %plot(f2,t(x),acc)
%     warning('off')
%     fitIac = struct(f2); fitIacv = struct(f2v);
%     warning('on')
%     Axc = fitIac.coeffValues{1,1};
%     Txc = fitIac.coeffValues{1,2};
%     Y0xc = fitIac.coeffValues{1,end};
%     R2xc = fitIacv.rsquare;
%     %fprintf('Ttotal(theory)=%1.2fms\n',Ttotaltheory)
%     Ttotalsimxc = Txc*1000;
    
    %fprintf('Ttotal(sim)=%1.2fms (R2=%1.2f)\n',Ttotalsim,R2ac)
    if photobleach==1
        paraout(end+1,:)=[Ttotalsimnac R2nac Ttotalsimacft R2acft Y0nac+Anac Y0acft+Aacft (Y0nac+Anac2)/Anac (Y0acft+Aacft2)/Aacft];
    elseif photobleach==0
        paraout(end+1,:)=[Ttotalsimnac R2nac Ttotalsimacft R2acft Y0nac+Anac Y0acft+Aacft (Y0nac)/Anac (Y0acft)/Aacft];
    end
end

fprintf('Input parameters:\n');
if photobleach ==1
    disp('photobleaching added')
end
fprintf('tau_on = %.4g ms, tau_off = %.4g ms, on Keq = %.4g \n',tau_on*1000,tau_off*1000,Keq_on);
fprintf('quenching efficiency = %.4g, photons/frame %.4g, offset %.4g \n',QE,emrate*dt,offset);
fprintf('predicted tau_tot = %1.2f, rough est. number of blinks %.4g\n\n',1/ktot*1000,nblink_est);

fprintf('fitting range %.4g to %.4g\n',t(min(ran)),t(max(ran)))

r2lim=0.9;

nacout=paraout(paraout(:,2)>r2lim,1);
fprintf('tau_tot(nac)av=%1.2f(%1.2f)ms, N(R2>0.9)=%2.0f\n',mean(nacout),std(nacout),sum((paraout(:,2)>r2lim)))
fprintf('nac(0)=%1.2e(%1.2e)\n',mean(paraout(:,5)),std(paraout(:,5)))
% [nac1,nac2]=hist(paraout(paraout(:,6)<5,6));
% naceq=fit(nac2(:),nac1(:),'gauss1');
% fprintf('k_eq=%1.2f\n',naceq(1))
% [nT, xoutT] = hist(nacout,10);
% [gT,gTv]=fit(xoutT(:), nT(:), 'gauss1');
% warning('off')
% fitgT = struct(gT); fitgTv = struct(gTv);
% warning('on')
% gMean = fitgT.coeffValues{1,2};
% gStd = fitgT.coeffValues{1,3}/sqrt(2);
% fprintf('tau_tot(nac)gauss=%1.2f(%1.2f)ms (R2>0.9)\n',gMean,gStd)
% [~,poiss]=poissfit(hist(paraoutfil(:,1)));
% fprintf('tau_tot(nac)poiss=%1.2fms (R2>0.9)\n',poiss(1))
%fprintf('tau_tot(nac)poly=%1.2f(%1.2f)ms (R2>0.9)\n',mean(paraout(:,7)),std(paraout(:,7)))


acftout=paraout(paraout(:,4)>r2lim,3);
fprintf('tau_tot(acft)av=%1.2f(%1.2f)ms N(R2>0.9)=%2.0f\n',mean(acftout),std(acftout),sum((paraout(:,4)>r2lim)))
fprintf('acft(0)=%1.2e(%1.2e)\n\n',mean(paraout(:,6)),std(paraout(:,6)))
%acftkeq=mean(paraout(paraout(:,6)<5,6));
%acfteq=fit(acfto2(:),acfto1(:),'gauss1');
%fprintf('k_eq=%1.2f N=%2i\n',acftkeq,sum(paraout(:,6)<5))
% [nT, xoutT] = hist(acftout,10);
% [gT,gTv]=fit(xoutT(:), nT(:), 'gauss1');
% warning('off')
% fitgT = struct(gT); fitgTv = struct(gTv);
% warning('on')
% gMean = fitgT.coeffValues{1,2};
% gStd = fitgT.coeffValues{1,3}/sqrt(2);
% fprintf('tau_tot(acft)gauss=%1.2f(%1.2f)ms (R2>0.9)\n',gMean,gStd)
% [~,poiss]=poissfit(hist(paraoutfil(:,1)));
% fprintf('tau_tot(acft)poiss=%1.2fms (R2>0.9)\n',poiss(1))

out(end+1,:)=[QE Keq_on mean(paraout(:,5))];

% xcout=paraout(paraout(:,6)>r2lim,5);
% fprintf('tau_tot(xc)av=%1.2f(%1.2f)ms N(R2>0.9)=%2.0f\n',mean(xcout),std(xcout),sum((paraout(:,4)>r2lim)))
% [nT, xoutT] = hist(xcout,10);
% [gT,gTv]=fit(xoutT(:), nT(:), 'gauss1');
% warning('off')
% fitgT = struct(gT); fitgTv = struct(gTv);
% warning('on')
% gMean = fitgT.coeffValues{1,2};
% gStd = fitgT.coeffValues{1,3}/sqrt(2);
% fprintf('tau_tot(xc)gauss=%1.2f(%1.2f)ms (R2>0.9)\n',gMean,gStd)
% [~,poiss]=poissfit(hist(paraoutfil(:,1)));
% fprintf('tau_tot(xc)poiss=%1.2fms (R2>0.9)\n',poiss(1))

% accun = acft(simcmos(1:20000));
% %accun = xcorr(simcmos(1:20000)-min(simcmos(1:20000)));accun=accun(20000+1:end);
% %accun = xcorr(simcmos(1:20000));accun=accun(20000+1:end);
% ran = 2:40;
% accun = accun(ran);
% accun = accun/max(accun);
% x = 1:length(accun);
% [fun,funv] = fit(t(x)',accun','A*exp(-x/T)+y0','startpoint',[1,0.01,1]);
% subplot(3,2,[5 6])
% plot(fun,t(x),accun)
% warning('off')
% fitIacun = struct(fun); fitIacunv = struct(funv);
% warning('on')
% Aun = fitIacun.coeffValues{1,1};
% Tun = fitIacun.coeffValues{1,2};
% Y0un = fitIacun.coeffValues{1,3};
% R2acun = fitIacunv.rsquare;
% %fprintf('QE(theory)=%1.2f\n',QE)
% fprintf('Keq(theory)=%1.2f\n',kon/koff)
% fprintf('Keq(sim)(AC)=%1.2f (R2=%1.2f)\n',Y0un/Aun,R2acun)
% fprintf('Keq(sim)(AC)=%1.2f (R2=%1.2f)\n',Tun/(Y0un/Aun),R2acun)

% %figure(2)
% [n1, xout1] = hist(simcmos,100);
% %baraxx=bar(xout1,n1); grid; hold off;
% [g1,g11]=fit(xout1(:), n1(:), 'gauss2','startpoint',[2000,150,10,2000,400,20]);
% warning('off')
% fitI = struct(g1); fitIv = struct(g11);
% warning('on')
% a1 = fitI.coeffValues{1,1};
% b1 = fitI.coeffValues{1,2};
% c1 = fitI.coeffValues{1,3};
% a2 = fitI.coeffValues{1,4};
% b2 = fitI.coeffValues{1,5};
% c2 = fitI.coeffValues{1,6};
% R2 = fitIv.rsquare;
% %fprintf('QE(sim)(H)=%1.2f\n',(b1-offset)/(b2-offset))
% fprintf('Keq(sim)(H)=%1.2f\n',(a2*c2)/(a1*c1))