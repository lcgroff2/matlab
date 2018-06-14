function results=etdiffnp_fun(parm)

% exdiff_fun -- excitons can diffuse, decay, or undergo energy transfer
% to an acceptor trap.  For simulating combination of energy
% diffusion and energy transfer inside a nanoparticle.

% Based on random walk simulation.  At each time step, the exciton
% can hop to a new site, decay, or undergo energy transfer.

% Usage: (See exdiffsetup.m for example)

% results = exhopfret_fun_opt2(parm)

% results.qyet = Energy Xfer efficiency (number of excitons
% that undergo energy transfer divided by the total number of excitons).

% results.decays = Fluorescence decays (list of lifetimes of
% individual excitons, can later make a histogram which can be compared to
% fluorescence decay data or other theories such as kww)

% Development history:

% based on exhopfret code, but this version does not use a a random walk
% on a cubic lattice, but rather, we set a certain timestep, and assume that
% the exciton changes position at each timestep .dt. by an amount:

% x = xprev + sign*randn(size(xprev));

% where

% sig = sqrt(2*D*dt);

% and ld = sqrt(2*D*tau), so D = ld^2/(2*tau)


%% End of simulation parameters

% all of the system physical parameters are in "SystemParams"
nprad = parm.nprad;
etrad = parm.etrad;
ld    = parm.ld;
tau0  = parm.tau0;
ndye  = parm.ndye;
xdye  = parm.xdye;
ydye  = parm.ydye;
zdye  = parm.zdye;

% All of the simulation-specific parameters are here.
% If chosen correctly, the results don't depend on these too much.
% This is where you sacrifice accuracy to gain speed, by choosing
% larger time steps, fewer excitons, or less averaging.
% So, you should try varying the values to make sure you haven't
% sacrificed too much accuracy.  You should make sure
% dt is less than 20*tau/(1-quenchingefficiency).  You should do
% several runs to check to make sure the results are consistent
% from run-to-run.  If not, then increase either nex or navg.
% Increasing nex requires more memory, so if more consistency is
% needed, and nex is already at say 3000, then increase navg.
nex  = parm.nex;
dt = parm.dt;
navg=parm.navg;

% diffusion parameters, based on ld and tau0
% diffusion constant
D = ld^2/(2*tau0);
% displacement per dt timestep
sig = sqrt(2*D*dt);

nsteps=round(tau0*6/dt); % do enough time steps so most excitons
                           % decay during the length of the random
                           % walk simulation
if nsteps<20,
  error('** WARNING ** number of timesteps too small, please decrease dt')
end


% Generate Random Exciton Positions Within Sphere of radius <nprad>
% pre-allocating arrays:
xpos=zeros([1 nex]);
ypos=xpos;
zpos=xpos;

cnt=1; % index counter

while cnt <= nex,
  xtry = nprad*2*(rand(1)-.5);
  ytry = nprad*2*(rand(1)-.5);
  ztry = nprad*2*(rand(1)-.5);
  if (xtry*xtry + ytry*ytry + ztry*ztry) < (nprad^2)
    xpos(cnt)=xtry;
    ypos(cnt)=ytry;
    zpos(cnt)=ztry;
    cnt=cnt+1;
  elseif parm.cubic == 1,
    xpos(cnt)=xtry;
    ypos(cnt)=ytry;
    zpos(cnt)=ztry;
    cnt=cnt+1;
  end
end


% arrays for keeping track of exciton fates
alive=zeros([1 nex])+1; % alive = 1 - decayed - transferred

transferred=zeros([1 nex]); % =1 for walkers that have
                                   % undergone FRET
decayed=zeros([1 nex]); % =1 for walkers that have decayed
                               % to ground state.
			       
decaytime=zeros([1 nex]); % record of time step at which a
                                 % given exciton decayed.

for curstep=1:nsteps,
  t=dt*curstep; % keep track of time
  % old: use random numbers to jump in either +/- X, Y, or Z direction
  
  %jumps = sign(rand([1 NumExcitons])-.5);
  
  % old: dir : 0 -> X hop, 1 -> Y hop
  
  %dir = floor(rand([1 NumExcitons])*3);
  
  % Some walkers decay to ground state
  justdecayed = alive .* (exp(-dt/tau0) < rand([1 nex]));
  alive = alive - justdecayed;
  decayed = decayed + justdecayed;

  % Record the decay time of each "justdecayed" exciton.
  decayidx = find(justdecayed==1);
  decaytime(decayidx) = t;
  
  % check to see if all excitons have decayed.  If so, exit loop
  nlive=sum(alive);
  aliveidx=find(alive);
  if nlive<1.
    break;
  end
  
  
  
  % Calculate overall energy transfer rate per timestep
  % (tricky!). Need to sum over walkers and dye acceptors.
  rates=zeros([1 nlive]);
  FR6T=etrad^6/tau0;
  for dyeidx=1:ndye,
    DX = xdye(dyeidx) - xpos(aliveidx);
    DY = ydye(dyeidx) - ypos(aliveidx);
    DZ = zdye(dyeidx) - zpos(aliveidx);
    Rsq = DX.*DX + DY.*DY + DZ.*DZ  + 0.0000001;
    rates = rates + FR6T ./(Rsq.*Rsq.*Rsq);
  end

  % Based on rate, use random lottery to determine which
  % walkers have undergone energy transfer.
  justtransferred = exp(-rates*dt) < rand([1 nlive]);
  alive(aliveidx)=alive(aliveidx)-justtransferred;
  transferred(aliveidx)=transferred(aliveidx)+justtransferred;
  
  % old code: only living walkers move
  %Xpos = Xpos + (jumps .* (dir==0)) .* alive*StepSize ;
  %Ypos = Ypos + (jumps .* (dir==1)) .* alive*StepSize ;
  %Zpos = Zpos + (jumps .* (dir==2)) .* alive*StepSize ;

  % new code, hops based on gaussian random number generator,
  % only 'living' excitons move.
  xhop = (sig*alive) .* randn([1 nex]);
  yhop = (sig*alive) .* randn([1 nex]);
  zhop = (sig*alive) .* randn([1 nex]);
  % add hops to previous positions
  xpos=xpos+xhop;
  ypos=ypos+yhop;
  zpos=zpos+zhop;
  
  % check for excitons that have escaped (R>ParticleRadius).
  % If there is an escape, then back off one step.
  
  if parm.cubic==0,
    escapelist=find((xpos.*xpos + ypos.*ypos + zpos.*zpos) > nprad^2);
  else
    if abs(xpos)>nprad, % this loop won't work as written, need to rewrite
      xpos=xpos-xhop;   % using 'find', if we want to acutally use cubic
    end
    if abs(ypos)>nprad,
      ypos=ypos-yhop;
    end
    if abs(zpos)>nprad,
      zpos=zpos-zhop;
    end
  end
  
  %if length(escapelist)>0,
  %  for eidx=escapelist,
  %    Xpos(eidx) = Xpos(eidx) - StepSize*sign(Xpos(eidx));
  %    Ypos(eidx) = Ypos(eidx) - StepSize*sign(Ypos(eidx));
  %    Zpos(eidx) = Zpos(eidx) - StepSize*sign(Zpos(eidx));
  %  end
  %end
  if length(escapelist)>0,
    for eidx=escapelist
      xpos(eidx) = xpos(eidx) - xhop(eidx);
      ypos(eidx) = ypos(eidx) - yhop(eidx);
      zpos(eidx) = zpos(eidx) - zhop(eidx);
    end
  end
  % double-check (comment out if it passes test)
  %escapelist=find((xpos.*xpos + ypos.*ypos + zpos.*zpos) > nprad^2);
  %if length(escapelist)>0
  %  disp('error in boundary-handling!!')
  %  error('boundary error!!');
  %end
end

% Calculate energy xfer efficiency as the fraction of walkers
% that have undergone fret:

results.qyet = (nex-sum(decayed)-sum(alive))/nex;

% record the decay times.
results.decays = decaytime(find(decayed));

% record some other stuff for debugging purposes
results.sig=sig;
results.D=D;
results.nlive=nlive;


% See if gives stretched-exponential behavior.

% Check to see if "local depletion" occurs correctly (few excitons
% in vicinity of the dyes).

% mean(results.decays) seems to give reasonable results (need to
% check with the detailed results, and see if reproduces
% stretch-exp behavior).

% Efficiency goes up if LD or ForsterRadius goes up.


