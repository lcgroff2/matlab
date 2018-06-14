
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
dt = 1e-3; % time per frame, in seconds

tmax=120; % number of seconds for experiment

tau_on = 0.055556; % lifetime of "on" state in seconds

tau_off = 0.5; % lifetime of "off" state in seconds

emrate=50000000; % number of photons per second, when in the "on" state

QE = 0.99; % quenching efficiency ( around 1.0 for single dye )


% camera settings and properties
offset=0; % zero if using some sort of baseline subtraction. for "raw" pixel data, use 100.
readnoise=2; % here we use units of adc counts, not photons, which is different than the way camera specs are given
gain=0.3; % gain of 0.3 means roughly 3 photons per ADC count

% simulation parameters (affect how things are calculated, but do not correspond to physical picture)
oversamp=10; % basically, uses a much smaller dt when calculating the trajectory, then samples back down to dt

% derived parameters, based on the above. Do not edit.
t=dt:dt:tmax; % time array

kon = 1 / tau_off; % kon = rate const for switching on

koff = 1/ tau_on ; % koff = rate const for switching off

Keq_on = kon / koff; % equilibrium constant, [on] / [off] or n_on / n_off


ktot = kon+koff; % total rate constant

% estimate of number of blink (on+off) events. This relationship is useful for designing power-dependent blinking experiments
nblink_est = tmax / (tau_on+tau_off); % see derivation in blinktstnotes.txt

% print out simulation parameters
fprintf(' input parameters:\n');
fprintf(' tau_on = %.4g ms, tau_off = %.4g ms, on Keq = %.4g \n',tau_on*1000,tau_off*1000,Keq_on);
fprintf(' quenching efficiency = %.4g, photons/frame %.4g \n',QE,emrate*dt);
fprintf(' predicted acorr tau = %.4g, rough est. number of blinks %.4g\n',1/ktot,nblink_est);


% use exponentially-distributed random numbers to generate on and off times
% should make a few extra, thus the factor of 1.5
ontimes = exprnd(tau_on,ceil(nblink_est*1.5),1);
offtimes = exprnd(tau_off, ceil(nblink_est*1.5),1);

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
    if curidx+nsteps>length(ts) % don't go past the end of "ts"
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

nblink=bidx-1;

fprintf('\nsimulation results:\n');
fprintf(' actual number of blink events in simulation: %i\n',nblink);

% resample molstate down so time steps are "dt"
tmp=0;
for idx=1:oversamp,
  tmp=tmp+osstate(idx:oversamp:end); % need to check this
end
molstate=tmp/10; % molstate is "1" for time steps where the molecule is "on", and "0" for timesteps where the molecule is "off", and some value in-between if the molecule switched during the timestep.


% Now generate simulated data

photonsperframe=emrate*dt;
% use approximation that for N photons average, replace with a gaussian random number with mu=N and sigma=N
% also take into account qe
onpart=mypoissrnd(QE*photonsperframe*molstate); % simulated additional photons when in "on" state
offpart=mypoissrnd((1-QE)*photonsperframe*ones(size(t))); % simulated photons when in "off" state and part of photons when in "on" state
simphotons=round(onpart+offpart); % round to an integer number of photons (should already be integer, but just in case I later change the above formulas).

% simulated cmos data, for one pixel. Includes offset, gain, and readnoise in simulation.
simcmos = offset+floor(gain*simphotons)+round(readnoise*(randn(size(t))));

% calculate a few correlation functions
acc=acft(simcmos); % unnormalized autocorrelation

nacc=nac(simcmos); % normalized G(tau)


% TODO xcorr(I^2,I)

% TODO print out some basic stuff
fprintf('acc(1) = %.5g, acc(2)=%.5g, acc(50) = %.5g, nacc(400)=%.5g\n',acc(1),acc(2),acc(50),acc(400));
%fprintf('acc ratio: %.5g\n',(acc(2)-acc(400))/acc(400));
fprintf('nacc(1) = %.5g, nacc(2)=%.5g, nacc(50) = %.5g, nacc(400)=%.5g\n',nacc(1),nacc(2),nacc(50),nacc(400));


%ran=1:100;
%est_tau=sum(nacc(ran).*t(ran))/sum(nacc(ran));
%fprintf('avg <tau> from acorr, from first 100 points: %.5g sec\n',est_tau)

%ran=1:200;
%est_tau=sum(nacc(ran).*t(ran))/sum(nacc(ran));
%fprintf('avg <tau> from acorr, from first 200 points: %.5g sec\n',est_tau);

ran=2:11; % don't fit the first point, since it can be mostly due to "white noise"
P=polyfit(t(ran),log(nacc(ran)),1); % semilog-plot analysis, only valid if nacc goes to zero at long times, and also only if nacc(ran) is above zero everywhere
est_tau=-1/P(1);
% second pass: if tau>>10*dt then we should include more data points in fit
if est_tau>6*10*dt
  ran=2:100;
  P=polyfit(t(ran),log(nacc(ran)),1); % semilog-plot analysis
  est_tau=-1/P(1);
end
fprintf('est tau from semilog fit to acorr: %.5g sec\n',est_tau);

% TODO do basic tau calculation using <tau> = sum(I*t)/sum(I), and also higher moments
% done: not really useful

% TODO do basic fitting (done).

