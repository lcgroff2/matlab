
% umbrella.m

% Uses Colbert/Miller's DVR method for solving the 1D Schrodinger equation,
% using Matlab.  It usually works well for things like
% the Harmonic Oscillator and harmonic double well, and the the
% Morse potential. For Coulombic potentials (H-like atom),
% a shooting method using adaptive Runge-Kutta, or similar, might work better
% with fewer evaluation points.

% Jason McNeill, January 2013 (based on dvrtest.m, February 2011)

% some needed constants
h= 6.62606876e-34; % J s
hbar = 1.0545716e-34 ; % J s
c = 2.99792458e8 ; % m/s
qe = 1.602176462e-19; % C
me = 9.1093819e-31; % electron mass in kg
na =  6.022142e+23; % avogadro, molecules per mole
mp=1.6726216e-27; % proton mass, kg

% *** start of user-defined parameters ***

% For umbrella-inversion potential of NH3
% calculation details
rmin=-1.2*1e-10; % Minimum distance (lower cutoff), in meters
rmax=1.2*1e-10; % maximum distance (upper cutoff), in Angstroms
npoints=1000; % number of points to use in calculation (memory use scales as n^2)
% potential details (from GAMESS, B3LYP-6-311G)
d=0.39839*1e-10; % displacement of N atom with respect to plane containing H atoms 
Vb=4.4644e-20; % inversion barrier height (in J)
hbar_wav=5.30884e-12; % hbar, in "wavenumber*seconds"
e1_wav=1131.77; % vibrational frequency (energy) at bottom of well (wavenumbers)
omega1=e1_wav/hbar_wav; % calculate angular frequency
e2_wav=892.32; % vibrational frequency (energy) at top of barrier (a negative value)
omega2=abs(e2_wav)/hbar_wav;
mu=1.18583; % reduced mass from normal mode analysis, relative to proton mass
m=mu*mp;
% give some points on curve, then fit a polynomial
dd=0.02*1e-10;
rp=[-dd, 0, dd, (d-dd), d, (d+dd)]; % first three are near barrier, last 3 near well
Vp= -0.5*m*omega2^2*rp.*rp+Vb;
Vp(4:6)=0.5*m*omega1^2*rp(1:3).*rp(1:3);
% fit a 5th-order polynomial to the curve
P=polyfit(rp,Vp,5);
r=linspace(rmin,rmax,npoints);
dr=abs(r(2)-r(1));
% now put the potential together piecewise
Vvec=polyval(P,r); % barrier piece: this part of the potential is valid near r=0, on the positive side
ridx=find(r>0.3*d); % the harmonic well piece
Vvec(ridx)=0.5*m*omega1^2*(r(ridx)-d).^2;
halfpoints=round(npoints/2); % slightly wrong if npoints is not even
Vvec(1:halfpoints)=Vvec(npoints:-1:(halfpoints+1)); % reflect about origin



% set up of DVR calculation
% Make a diagonal matrix with Vvec along the diagonal:
Vmat=diag(Vvec);

% we decompose the kinetic energy expression into 4 parts:
% Kconst = constant prefactor = hbar^2/(2*m*dr^2)
% Kdiag = pi^2/3 - 1 / (2*i^2) = the diagonal part of the KE
% Ksignmat = (-1)^(i-j) = (-1).^(XX-YY)
% Knondiag = 2/(XX-YY)^2 - 2/(XX+YY)^2 = the non-diagonal part of the KE

Kconst=hbar^2/(2*m*dr^2); % factor on KE term
xtmp=1:npoints;
[XX,YY]=meshgrid(xtmp,xtmp);
Ksignmat=(-1).^(XX-YY);
Kdiag=pi^2/3 - 0.5*xtmp.^(-2);
Knondiag = 2*(XX-YY).^(-2) - 2*(XX+YY).^(-2);
Kmat=Knondiag;
for x=xtmp, % replace the diagonals with the correct expr
  Kmat(x,x)=Kdiag(x);
end
Kmat=Kmat.*Ksignmat; % multiply the signs term
Kmat=Kmat*Kconst; % multiply by the constant
  
Hmat=Kmat+Vmat; % The Hamiltonian Matrix is the sum of the kinetic
                % energy matrix and the potential matrix

[wavefunmat,eigvalmat]=eig(Hmat);

energies=xtmp;
for x=xtmp,
  energies(x)=eigvalmat(x,x);
end

fprintf(1,'First 4 energies are ');
fprintf(1,'%5g, ',energies(1:4));
fprintf(1,'J, or\n');
fprintf(1,'%5g, ',energies(1:4)/h/c/100); 
fprintf(1,' in wavenumbers\nPlotting first 4 wavefunctions and potential w/ energies...\n');

figure(1)
plot(r*1e10,wavefunmat(:,1),'b');
hold on
plot(r*1e10,wavefunmat(:,2)+.1,'g');
plot(r*1e10,wavefunmat(:,3)+.2,'r');
plot(r*1e10,wavefunmat(:,4)+.3,'c');
hold off
xlabel('r, Angstroms')
ylabel('probability amplitude');
axis tight

figure(2)
plot(r*1e10,Vvec);
hold on
for x=1:6,
  ran=find(Vvec<energies(x));
  plot(r(ran)*1e10,ones(size(Vvec(ran)))*energies(x),'r');
end
hold off
xlabel('r, Angstroms')
ylabel('Energy, J')
axis tight 
