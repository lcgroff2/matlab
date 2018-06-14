
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

tau_on = 0.030; % lifetime of "on" state

tau_off = 0.030; % lifetime of "off" state

emrate=200000; % number of photons per second, when in the "on" state

QE = 0.05; % quenching efficiency ( around 1.0 for single dye )


% camera settings and properties
offset=100;
readnoise=3;
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

simcmos = offset+floor(gain*simphotons)+round(readnoise*(randn(size(t))));
figure(1)
subplot(2,2,[1 2])
plot(simcmos(1:1000))
acc = acft(simcmos(1:20000));
ran = 2:200;
acc = acc(ran);
acc = acc/max(acc);
x = 1:length(acc);
g = fittype('0.0001*a*exp(-x/(1000*b))+c');
f2 = fit(x',acc',g);
subplot(2,2,[3 4])
plot(f2,x,acc)
