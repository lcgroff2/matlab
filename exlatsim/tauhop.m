%Determining the rate/time constant for ET per picosecond for polymer
%CPNs.

%Define relevant parameters for hopping probability:

dt = 1; %time step set to 1 ps
tau = 300; %300 ps lifetime for MEH-PPV in good solvent
ld = 12; %12 nm exciton diffusion length
dx = 0.9; %chromophore size determined for MEH-PPV from latspace.m script
          %(2 rep units).
D = ld^2/2/tau; %exciton diffusion constant in nm^2/ps
phop = 2*D*dt/dx^2; %Exciton hopping probability for random walk

%Determine hopping rate constant from phop = 1-khop*dt
khop = (1-phop)/dt; %Rate constant units of ps^-1 
tauhop = 1/khop; %Should give rough estimate of correlation time from 
                 %anisotropy, units of ps.
fprintf(1,'Exciton hopping probability = %1.3g\n Hopping Rate Constant = %1.3e ps^-1\n Time per Excition Hop = %1.2g ps',phop,khop,tauhop) 