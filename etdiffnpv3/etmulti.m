
% script for doing multiple etdiffnp3 runs sequentially and saving
% each result, for a given number of dyes, in a separate matfile.
% The most important final results are in qevec (quenching efficiency 
% for each given number of dyes), taukwwvec, betavec.

% number of dyes to do simulations for.
dyevec=[1:10];

%example run:
% nohup octave -q etmulti.m > etmulti1.log 2>&1 &

% nprad = nanoparticle radius, nm:
parm.nprad = 4; 
% etrad = Forster energy transfer radius, nm:
parm.etrad = 3.5;
% ld = exciton diffusion length, nm:
parm.ld = 6;
% tau0 = lifetime, in absence of quencher, picoseconds:
parm.tau0 = 3000;
% ndye = number of dyes/quenchers (per particle):
parm.ndye = 1;
% set cubic=1 for a 'cubic' particle (not yet fully implemented)
parm.cubic = 0;
%parm.nex = 6000; % 500 is probably the minimum number to use.  3000 typical.
parm.nex=3000;
% dt = time step, in picoseconds.  Never use larger than tau0/20.
% accuracy is compromised for dt>tau0/100.
parm.dt = 1;
% navg = number of averages - at least 20 needed (guess), 500 typical
%parm.navg = 1000;  
parm.navg=50;

% fill up plist
for idx=1:length(dyevec)
  plist{idx}=parm;
  plist{idx}.ndye=dyevec(idx);
end

% storage for key results
betavec=zeros(size(dyevec));
taukwwvec=betavec;
tauavgvec=betavec;
qevec=betavec;
tau1vec=betavec;
tau2vec=betavec;
fvec=betavec;
tauvec=betavec;

  
dtstr=datestr(date,'yyyy-mm-dd');

for pidx=1:length(dyevec)
    etdiffnp3
    fitetdiff
    savnam=sprintf('et-%s-%02i.mat',dtstr,pidx);
    clear ans myfun
    save('-mat',savnam);
end

