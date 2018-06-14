%Determining the lower limit for the fill fraction for percolation model
%calculations.

%Start from PS data for swelling. We extrapolated from a linear fit of the
%PS data that for 100% THF, there would be ~77% increase in lattice
%spacing.
sizeincrease = 1.77;

cubediam = 10; % 10 nm cube
cubevol = cubediam^3; % 10^3 nm^3
cubediam2 = cubediam*sizeincrease; %max particle size increase from swelling sims, nm.
cubevol2 = cubediam2^3; %nm^3

vacancy = cubevol2-cubevol; %difference in particle volumes should equate 
%to the vacancy space in the nanoparticle, assuming that the unswelled
%nanoparticle has no empty space between chromophores.

minfraction = 1-vacancy/cubevol2; %the amount of filled space is 1-(vacant volume/total swelled particle volume).

fprintf(1,'Lower limit on CPN fill fraction = %1.2g\n',minfraction);