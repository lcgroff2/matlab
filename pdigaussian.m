
% PDIGAUSSIAN.M - Determine the polydispersity index, based on a
% gaussian size distribution. There is probably an analytical solution,
% but I couldn't find it quickly so it seemed easier to make a script
% using a for-loop to estimate the PDI if we know the mean and std deviation
% in the particle size, from imaging data (AFM). Since this is from a particle
% histogram, it is assumed that the inputs are "number-averaged".

rmean = 8 ; % mean size (diameter or radius), from histogram

rstd = 2 ; % standard deviation (sigma), from histogram

dr = rstd/20;

r=0:dr:(rmean+4*rstd);

numdensity=(1/rstd/sqrt(2*pi))*exp(-(r-rmean).^2/2/rstd^2);

% weight-density is the number density, weighted by the MW (heavier particles count more)

wtdensity = numdensity.*r.^3;

% normalize wtdensity to 1 as needed for any PDF

wtdensity=wtdensity/sum(wtdensity*dr);

% find weight-averaged size

r_wt = sum(wtdensity.*r*dr);

pdi = r_wt / rmean;

fprintf(1,'for histogram yielding average size of %4g and std dev of %4g\n', rmean, rstd);
fprintf(1,'the mass-weighted average size is %4g and the PDI is %4g\n',r_wt, pdi);

