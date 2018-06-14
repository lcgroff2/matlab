

% etdiffnp3.m

% script for simulating enrgy transfer and exciton diffusion in a nanoparticle.

% the main script is "etmulti".  This is a helper script, 

% make a matrix for storing ET efficiency results

if exist('pidx')~=1,
  disp('run etmulti instead')
  return
end

parm=plist{pidx};
parm
rand('state',sum(100*clock)); % randomize the random number generator



% Do some sanity-checking here
D = parm.ld^2/(2*parm.tau0);
% displacement per dt timestep
sig = sqrt(2*D*parm.dt);
if sig>parm.etrad/10,
  disp('warning: rms displacement per timestep > parm.etrad/10 , consider decreasing parm.dt');
end
if parm.navg<40,
  disp('warning, possibly not enough averaging, consider increasing parm.navg');
end
if parm.dt>parm.tau0/50,
  disp('warning: parm.dt>parm.tau0/50, consider increasing parm.dt');
end


effmat=zeros([1 parm.navg]);

decaylist=[];

for avcnt=1:parm.navg,  % average over many particles (many sets of
                        % dye positions)
  ndye=parm.ndye;
  nprad=parm.nprad;
  xdye=zeros([1 parm.ndye]); % pre-allocate empty arrays
  ydye=xdye;
  zdye=xdye;


  % loop to generate random dye positions within the particle radius
  cnt=1; % index counter
  while cnt <= parm.ndye,
    xtry = parm.nprad*2*(rand(1)-.5);
    ytry = parm.nprad*2*(rand(1)-.5);
    ztry = parm.nprad*2*(rand(1)-.5);
    rsq = xtry*xtry + ytry*ytry + ztry*ztry;
    if (rsq < parm.nprad^2)
      xdye(cnt)=xtry;
      ydye(cnt)=ytry;
      zdye(cnt)=ztry;
      cnt=cnt+1;
    elseif parm.cubic == 1,
      xdye(cnt)=xtry;
      ydye(cnt)=ytry;
      zdye(cnt)=ztry;
      cnt=cnt+1; 
    end       
  end
  meanradius=mean(sqrt(xdye.*xdye + ydye.*ydye + zdye.*zdye));
  % store a copy of xdye, ydye, zdye
  parm.xdye=xdye;
  parm.ydye=ydye;
  parm.zdye=zdye;
  % Do a simulation for one particle
  %return; % for testing
  results = etdiffnp_fun3(parm);
  % Store quantum yield of energy transfer data in matrix "qyet".
  qyet(avcnt) = results.qyet;

  % Store the list of all the decay times.  This is related to
  % TCSPC results.
  % decaylist=cat(2,decaylist,results.decays(:)');
  if avcnt==1,
      fdecay=results.population;
  else
      fdecay=fdecay+results.population;
  end
end




%tt=histran;
%fdecay=hist(decaylist,tt);
%fdecay(end)=fdecay(end-1); % remove spike at the end due to way hist function works
%plot(tt,fdecay)
tt=results.time;
xlabel('decay time');
ylabel('number of excitons');
title('fluorescence lifetime simulation (decay time histogram)');

fprintf(1,'Quenching efficiency (1.0=max): %6.4f\n',mean(qyet));
qe=mean(qyet);
fprintf(1,'Uncertainty in quenching eff: %8.5f\n',std(qyet)/sqrt(parm.navg));
fprintf(1,'RMS displacement per time step: %8.5f nm\n',results.sig);

pt=fdecay/sum(fdecay);
halflife=sum(pt.*tt);
fprintf(1,'Mean exciton lifetime %7.1f ps\n',halflife/log(2));

if results.sig>(parm.etrad/5),
  fprintf(1,'warning: sig too large compared to Forster radius, consider decreasing dt\n');
end

if parm.dt>halflife/40,
  fprintf(1,'Warning: dt > lifetime/40. Consider decreasing dt\n')
end
if (mean(qyet)>0.05) && (std(qyet)/mean(qyet)/sqrt(parm.navg)>0.1),
  fprintf(1,'Warning: relative uncertainty seems high, consider increasing navg\n')
end

%** beta estimation, using 'log-log-log' method **
% see "logloglog.pdf" for details.
lt=log(tt);
lf=-log(fdecay/max(fdecay));
%[mx,idx]=max(fdecay);
idx=find(fdecay==max(fdecay));
if length(idx)>1, % check for the case of 2 or more points=max
  idx=max(idx);
end
lt=lt((idx+1):end); % remove some early time points
lf=lf((idx+1):end);
% remove late time points where fdecay=0;
llf=log(lf);
idx=find(isinf(llf));
if length(idx)>0,
  minidx=min(idx);
  lt=lt(1:(minidx-1));
  llf=llf(1:(minidx-1));
end
pp=polyfit(lt,llf,1);

fprintf(1,'estimated beta (approx) %7.4f\n',pp(1));

%fprintf(1,'type "parm" to see/check input parameters\n');
fprintf(1,'fluorescence lifetime simulation in "tt" and "fdecay", plot(tt,fdecay)\n');

%qevec(pidx)=mean(qyet);

