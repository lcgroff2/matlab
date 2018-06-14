
% exlatsim.m

% lattice gas simulation, of excitons in a "swelled" (swollen?) nanoparticle or polymer chain in solution

% this version assumes a hopping rate based on the exciton diffusion constant

% JDM March 2014

% simulation parameters -- in later versions these will be specified in a control script instead
dx = 0.9 ; % lattice spacing, nm
ld = 12 ; % exciton diffusion length, nm
tau = 320; % lifetime in solution, ps
dt = 1; % time step, ps
dxsw = 1.60; % Swelled lattice spacing, nm
k = 4; % hopping probability depends on swelled spacing to an exponent k (between 2-4)
%diaghop=0; % allow diagonal "hop" motion? (probability should be less). NOT IMPLEMENTED YET
cubesz=round(10*(dxsw/dx)); % 10 nm cube initially, diameter (and thus, the lattice) increases as diameter*f^(-1/3)
fillfrac=1.0; % chromophore filling fraction (1.0 = all sites filled)
R0=4; % forster radius
numq=10; % number of quenchers 
numex=1000; % number of excitons

% derived parameters
D1D = ld^2/6/tau; % 1D exciton diffusion constant 
%phop=0.1; % need to replace with proper expression
% for random-walk with a hopping probability, sigma^2 = p (t/dt) dx^2 = 2 D t (for 1D),
% therefore phop = 2 D dt / dx^2. Accounting for lattice swelling, multiply
% by swelling fraction f^(k/3)
phop=2*D1D*dt/dx/dx;
% %Uncomment below for exponentially decaying phop.
% kappa = log(1/phop)/dx;
% phop = exp(-kappa*dxsw);
% Uncomment below to adjust phop by swelling fraction.
fk3 = (dx/dxsw)^k; %Swelling fraction f^(k/3), positive exponent because otherwise phop increases as dxsw increases.
phop = phop*fk3;
fprintf(1,'swelling fraction (0-1): %1.3g\n',fk3);
% phopdiag=phop/2; % by geometry, diagonal hops are less probable (not implemented)
fprintf(1,'hopping probability: %.3g per timestep\n',phop);
fprintf(1,'make sure it is neither too small nor too large: 1-10%% is probably OK\n');

% make new random seed using system clock
rand('state',sum(100*clock));

% make cubic matrix representing which sites are occupied by chromophores
gridpoints=round(cubesz/dxsw); % number of grid points along each coordinate axis

% occupy sites at random. But in future could use a self-avoiding random walk, etc.
cmat=rand([gridpoints gridpoints gridpoints])<fillfrac;

% make XX,YY,ZZ meshgrid for cubic particle
x=(1:gridpoints)*dx;
[XX,YY,ZZ]=meshgrid(x,x,x);

% make quenchers--for now make them random numbers not on site and not on chromophore. It shouldn't
% make much difference for the random lattice gas model. But for other systems (random coil, etc)
% we will need to change it.
if numq>0
  qposx=dx+rand([1 numq])*(gridpoints-1)*dxsw; % need to be careful of off-by-one errors here, since we need to make sure the quenchers overlap with where the chromophores are.
  qposy=dx+rand([1 numq])*(gridpoints-1)*dxsw; 
  qposz=dx+rand([1 numq])*(gridpoints-1)*dxsw; 
end

% % for testing: put a single quencher at dead center
% qposx=cubesz/2; qposy=qposx; qposz=qposx; disp('*** TESTING w/ 1 QUENCHER AT DEAD CENTER ***');

% exciton indices, exi, exj, exk, random number from 1 to gridpoints
exi=zeros([1 numex]);
exj=exi;
exk=exi;
for idx=1:numex,
  goodpos=0;
  while goodpos==0
    tryi=1+floor(rand(1)*gridpoints);
    tryj=1+floor(rand(1)*gridpoints);
    tryk=1+floor(rand(1)*gridpoints);
    if cmat(tryi,tryj,tryk), % check to see if chromophore exists at that point
      exi(idx)=tryi;
      exj(idx)=tryj;
      exk(idx)=tryk;
      goodpos=1;
    end
  end
end        

% record initial positions (needed for debugging, and twinkling simulations)
exini=exi;
exinj=exj;
exink=exk;


% prepare for main loop
nsteps=round(6*tau/dt);
expop=zeros([1 nsteps]); % store exciton population here
exstate=ones([1 numex]); % exciton state=1 for alive exciton, =0 for decayed
t=(1:nsteps)*dt;
pdec=1-exp(-(dt/tau)); %probability of decay by rad/nonrad transition during 1 timestep
net=0; % number of excitons that have undergone energy xfer
ndec=0; % number of excitons that have decayed by rad or nrad relaxation
etidxstore=[]; % store the indices of the excitons that undergo ET
decidxstore=[]; % store the indices of the excitons that undergo regular rad/nonrad decay

for tidx=1:nsteps, % loop over time steps
  % first find which have decayed by fluorescence and non-radiative recomb,
  % also taking into account which excitons have already decayed.
  exdecayed=exstate.*(rand([1 numex])<pdec);
  % update state (1=alive, 0=dead)
  exstate=exstate-exdecayed;
  decidxstore=[decidxstore find(exdecayed)]; % store indices of excitons that decayed
  ndec=ndec+sum(exdecayed); % update number of excitons that have undergone decay
    
  % calculate ET rate for each exciton (to each quencher)
  ket=zeros([1 numex]); % zero out because some excitons have decayed
  if numq>0 % only do energy xfer calculation if there are quenchers/dyes
    for eidx=1:numex % TODO: in optimized version, will only loop over live excitons and eliminate "if" statement below
      ktot=0; % total energy transfer rate for this exciton
      if exstate(eidx) % only do calculation if exciton still alive, to save computation time
        %for qidx=1:numq % might try to vectorize this, later
        %Rsq=(dx*exi(eidx)-qposx(qidx))^2 + (dx*exj(eidx)-qposy(qidx))^2 + (dx*exk(eidx)-qposz(qidx))^2;
        %ktot=ktot+(1/tau)*R0^6/Rsq/Rsq/Rsq; % decay rate from FRET theory
        %end
        % here is vectorized version
        Rsq=(dxsw*exi(eidx)-qposx).^2 + (dxsw*exj(eidx)-qposy).^2 + (dxsw*exk(eidx)-qposz).^2;
        ket(eidx)=sum((1/tau)*R0^6./(Rsq.*Rsq.*Rsq));
        %ket(eidx)=ktot; % the total energy transfer rate for this exciton
      end
    end
  
    % now determine which excitons decay via energy transfer
    pet=1-exp(-dt*ket); % probability of ET
    exdecayed=rand([1 numex])<pet; % don't need to multiply by exstate, because pet is zero for those
    exstate=exstate-exdecayed;
    etidxstore=[etidxstore find(exdecayed)]; % store indices of excitons that did ET
    net=net+sum(exdecayed); % update number of excitons that have undergone ET
  end
  
  % update population
  expop(tidx)=sum(exstate);
  
  
  % now allow excitons to hop, but only if there is a chromophore there and it is live
  % first for x axis
  chkstep = rand([1 numex])<phop;
  chkstep = chkstep.*exstate; % check if alive
  
  trymove=floor(rand([1 numex])*2)*2-1; % this gives +1, -1 w/ 50-50 probability
  trymove=trymove.*chkstep; % only alive excitons with enough probability can move
  
  % check to see if there is a chromophore there (should try to vectorize this)
  for midx=find(trymove ~=0), % only iterate over actual moves   
    % test to make sure not out of bounds
    if ((trymove(midx)+exi(midx))>0) && ((trymove(midx)+exi(midx))<gridpoints+1)
      if cmat(trymove(midx)+exi(midx),exj(midx),exk(midx))==1,
        exi(midx)=exi(midx)+trymove(midx);
      end
    end
  end
  
  % next, y axis
  chkstep = rand([1 numex])<phop;
  chkstep = chkstep.*exstate; % check if alive
  
  trymove=floor(rand([1 numex])*2)*2-1; % this gives +1, -1 w/ 50-50 probability
  trymove=trymove.*chkstep; % only alive excitons with enough probability can move
  
  % check to see if there is a chromophore there (should try to vectorize this)
  for midx=find(trymove ~=0), % only iterate over actual moves
    % test to make sure not out of bounds
    if ((trymove(midx)+exj(midx))>0) && ((trymove(midx)+exj(midx))<gridpoints+1)
      if cmat(exi(midx),trymove(midx)+exj(midx),exk(midx))==1,
        exj(midx)=exj(midx)+trymove(midx);
      end
    end
  end
  
  % finally, z axis
  chkstep = rand([1 numex])<phop;
  chkstep = chkstep.*exstate; % check if alive
  
  trymove=floor(rand([1 numex])*2)*2-1; % this gives +1, -1 w/ 50-50 probability
  trymove=trymove.*chkstep; % only alive excitons with enough probability can move
  
  % check to see if there is a chromophore there (should try to vectorize this)
  for midx=find(trymove ~=0), % only iterate over actual moves
    if ((trymove(midx)+exk(midx))>0) && ((trymove(midx)+exk(midx))<gridpoints+1)
      if cmat(exi(midx),exj(midx),trymove(midx)+exk(midx))==1,
        exk(midx)=exk(midx)+trymove(midx);
      end
    end
  end
  
end

qe=net/numex;

%Estimate beta:
fdecay = expop/sum(expop);
lt = log(t);
lf = -log(expop/max(expop));
idx=find(fdecay==max(fdecay));
if length(idx)>1, % check for the case of 2 or more points=max
    idx=max(idx);
end
lt=lt((idx+1):end); % remove some early time points
lf=lf((idx+1):end);
llf = log(lf);
idx=find(isinf(llf));
if length(idx)>0,
    minidx=min(idx);
    lt=lt(1:(minidx-1));
    llf=llf(1:(minidx-1));
end
pp=polyfit(lt,llf,1);
fprintf(1,'estimated beta (approx) %7.4f\n',pp(1));

tauave=sum(t.*expop)/sum(expop);
fprintf('System parameters: %i quenchers, LD=%3g nm, R0=%3g nm, tau=%3g ps, cube size=%3g nm, fill=%4g\n', numq, ld, R0, tau, cubesz,fillfrac);
fprintf('sim parameters: dt=%3g ps, dx=%3g nm, dxsw=%3g nm\n',dt,dx,dxsw);
fprintf('predicted occupied volume: %.4g nm^3\n',fillfrac*cubesz^3);
fprintf('actual occupied volume: %.4g nm^3\n',sum(sum(sum(cmat)))*dxsw^3);
fprintf('average/expected number of occupied sites: %i\n',round(fillfrac*gridpoints^3));
fprintf('actual number of of occupied sites: %i\n',sum(sum(sum(cmat))));
fprintf('Results: QE=%4g, tauave=%4g\n',qe,tauave);

% plot(t,expop)

% make list of exciton positions for later analysis (convert etidxstore, decidxstore into positions)
% this is used by pltrdist
if numq >= 1,
etposstore=zeros([3 length(etidxstore)]);

etposstore(1,:)=exini(etidxstore)*dxsw;
etposstore(2,:)=exinj(etidxstore)*dxsw;
etposstore(3,:)=exink(etidxstore)*dxsw;
end

decposstore=zeros([3 length(decidxstore)]);

decposstore(1,:)=exini(decidxstore)*dxsw;
decposstore(2,:)=exinj(decidxstore)*dxsw;
decposstore(3,:)=exink(decidxstore)*dxsw;

