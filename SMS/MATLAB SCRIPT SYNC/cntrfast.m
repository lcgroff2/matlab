function centers = cntrfast(img,pks,sz)

% finding peak positions using custom algorithm for calculating the
% centroid quickly.

% **NOTE** for this fit routine, we must run bpass or fpass on image first

% **also note that the Grier/Crocker routine uses a circular mask, while 
% we do not, but it should be easy enough to add, just by zeroing out
% any pixels outside the circle

% **also note that this function only works for a single "pks" peak, not multiple
% peaks as cntrd does

% **at some point we redefined "sz" so it is equivalent to "sz/2" in Grier/Crocker routines

% **note that passing "img" to this function might take many cpu cycles, so should consider
% folding this function into the tracking script, for speed, or modifying so that many peaks
% are analyzed simultaneously, or making "img" global.

% Set up X and Y variables
% Be careful not to mix up rows/columns (need to check later)
% imgsz=size(img);
% X=1:imgsz(2); % let X be rows (need to check later)
% Y=1:imgsz(1); % let Y be columns

x0=pks(1); % center of gaussian from "peaks"
y0=pks(2);


% use only points within ROI for later calculations
XROI=round(x0-sz):round(x0+sz);
YROI=round(y0-sz):round(y0+sz);
imgROI=img(YROI,XROI);


% make vectors for x, y
% 
% x=XROI(:);
% y=YROI(:);

% sum intensities over rows and columns within Region Of Interest
IXROI=sum(imgROI);
%IXROI=IXROI(:); % make it a column vector
IYROI=sum(imgROI'); % rotate imgROI by 90' and then sum

%IYROI=IYROI(:);

% calculate <x>, <y>, <x^2>, <y^2>
Itot=sum(IXROI);
xbar=sum(IXROI.*XROI)/Itot;
ybar=sum(IYROI.*YROI)/Itot;
xsqbar=sum(IXROI.*XROI.*XROI)/Itot;
ysqbar=sum(IYROI.*YROI.*YROI)/Itot;

% calculate variance, sigma^2 (sometimes called radius of gyration squared in this context)
varx=xsqbar-xbar^2;
vary=ysqbar-ybar^2;
varave=(varx+vary)/2;

centers=[xbar ybar varave Itot];