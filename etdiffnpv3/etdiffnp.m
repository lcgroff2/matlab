
% etdiffnp.m -- "energy transfer and exciton diffusion in a nanoparticle"

% 3/2012 JDM

% script for calculating quenching efficiency and fluorescence
% dynamics, in nanoparticles doped with quenchers, such as dyes or
% polarons/defects.

% uses/requires "etdiffnp_hlp.m" and "etdiffnp_fun.m".

% This script is primarily for setting the simulation parameters
% (below).  "etdiffnp_hlp" does some additional setup (such as
% placing dyes in random positions), calls "etdiffnp_fun" (which
% does the actual simulation) and does some analysis and
% plotting of results returned by "etdiffnp_fun".  The actual
% random walk trajectories are generated in "etdiffnp_fun", along
% with all of the physics.  This script and "etdiff_hlp" only do
% setup and analysis.

% main simulation parameters, EDIT THESE TO MATCH YOUR SYSTEM.

% nprad = nanoparticle radius, nm:
parm.nprad = 12; 
% etrad = Forster energy transfer radius, nm:
parm.etrad = 3;
% ld = exciton diffusion length, nm:
parm.ld = 12.6;
% tau0 = lifetime, in absence of quencher, picoseconds:
parm.tau0 = 3000;
% ndye = number of dyes/quenchers (per particle):
parm.ndye = 22;
% set cubic=1 for a 'cubic' particle (not yet fully implemented)
parm.cubic=0;

% some notes about the parameters below:
% The number of excitons shouldn't affect the underlying physics,
% but does affect the accuracy of the result--you need a lot of excitons
% to properly sample the volume of the particle and to beat the
% poisson noise inherent in Monte Carlo simulations.
% Regarding the other parameters below here:
% If chosen correctly, the results don't depend on these too much.
% This is where you sacrifice accuracy to gain speed or vice-versa,
% by choosing larger time steps, fewer excitons, or less averaging.
% So, you should try varying the values to make sure you haven't
% sacrificed too much accuracy.  You should make sure:

% dt < 20*tau0/(1-quenchingefficiency)

% In other words, if there is a lot of quenching, then the average
% exciton lifetime is significantly less than tau0, so when this is
% the case, you should decrease dt so there are enough time steps
% to get a meaningful fluorescence lifetime simulation.

% You should do several runs to check to make sure the results are
% consistent from run-to-run.  If not, then increase either nex or
% navg.  Increasing nex requires more memory, so if more consistency
% is needed, and nex is already at say 3000, then increase navg.
% navg should be at least 20-50, to properly sample the space of 
% possible dye positions, and nex should probably be at least 500, to
% properly sample all of the possible exciton positions.

% nex = number of excitons:
parm.nex = 3000; % 500 is probably the minimum number to use.  3000 typical.
% dt = time step, in picoseconds.  Never use larger than tau0/20.
% accuracy is compromised for dt>tau0/100.
parm.dt = 1;
% navg = number of averages - at least 20 needed (guess)
parm.navg = 50;  
% histogram binning range for plotting fluorescence lifetime simulation.
histran=10:2:9000; % 10:20:9000 means histogram bins spaced 20 ps apart.

%% end of adjustable parameters.

% override above parameters if plist is given
if exist('plist','var')
  parm=plist{pidx};
end

% output parameters
parm

% give random number generator a seed based on current time, to ensure
% that we get a different set of random numbers each time we run.
rand('state',sum(100*clock));

tic;
etdiffnp_hlp
toc
% edit the line below if you want to save to a different filename
save etdiffnpcalc.mat

if exist('plist','var')
  qevec(pidx)=mean(qyet);
end

% little hack to work with "mlspawn.sh" script (ignore):
%if exist('instance')==1,
%  savnam=sprintf('exdiffresults%02i',instance);
%else
%  savnam='exdiffresults';
%end