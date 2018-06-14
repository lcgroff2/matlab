%hoprate.m, LCG 7/24/2014
%For determination the rate/time constant for ET per picosecond for polymer CPNs.

%Define relevant parameters for hopping probability:

dt = 1; %time step set to 1 ps
tau = 290; %300 ps lifetime for MEH-PPV in good solvent
ld = 12; %12 nm exciton diffusion length
dx = 0.9; %chromophore size determined for MEH-PPV from latspace.m script
          %(2 rep units).
D = ld^2/6/tau; %1D exciton diffusion constant in nm^2/ps
phop = 2*D*dt/dx/dx; %Exciton hopping probability for random walk

%Determine hopping rate constant from phop = 1-exp(-khop*dt), phop ~ khop/dt
khop = phop/dt; %Rate constant units of ps^-1 
tauhop = 1/khop; %Should give rough estimate of correlation time from 
                 %anisotropy, units of ps.
fprintf(1,'Exciton hopping probability = %1.3g\nHopping Rate Constant = %1.2g ps^-1\nTime per Excition Hop = %1.2g ps\n\n',phop,khop,tauhop) 