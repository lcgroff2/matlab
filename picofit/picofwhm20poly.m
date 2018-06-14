% function for finding FWHM of IRF for pico data
% some adjustable parameters
rep_rate=76.082e6;
period=1/rep_rate*1e9;
bins_per_ns = 377.97;
numpts=100; % number of points to include in centroid calc.

cd c:\matlab\picofit
ls *.asc
% [fnam,fpath]=uigetfile('*.asc','MCA Data');

fnam=input('enter file name ','s');
fid=fopen(fnam,'r');
y=fscanf(fid,'%f',[1,inf]);
fclose(fid);

x=1:length(y);
plot(x,y)
disp('click on a peak to measure FWHM');
plot(x,y)
[xin,yin]=ginput(1);
tmpmax=round(xin);

%Find max close to cursor
srchran=80;
[val,mxidx]=max(y(tmpmax-srchran:tmpmax+srchran));
%mxidx=mxidx+tmpmax-srchran;   
tmpx=x(tmpmax-srchran:tmpmax+srchran);
mxidx=tmpx(mxidx);
%Fit a parabola over the region close to the peak
ran=mxidx-2:mxidx+2;
p=polyfit(x(ran),y(ran),3);
a=p(1); b=p(2); c=p(3); d=p(4);
xmax=(-2*b)+sqrt(((2*b)^2)-(12*a*c))/(6*a);
ymax=a*xmax^3+b*xmax^2+c*xmax+d;
%Find point closest to the true left and right half-heights
lran=mxidx-srchran:mxidx;
uran=mxidx:mxidx+srchran;
lfwhm=min(find(y(lran)>ymax/2))+mxidx-srchran;	%This will give the nearest point above the true half max
ufwhm=min(find(y(uran)>ymax/2))+mxidx;
%Select the three points above and below the true max
lfwhmran=lfwhm-3:lfwhm+2;
ufwhmran=ufwhm-2:ufwhm+3;
%Fit straight lines to the parts of the curve near the left and right true half max locations
lp=polyfit(x(lfwhmran),y(lfwhmran),1);
up=polyfit(x(ufwhmran),y(ufwhmran),1);
%Determine where the fits intersect the true half max and determine the fwhm
fwhm=((ymax/2)-up(2))/up(1)-((ymax/2)-lp(2))/lp(1);

fwhm_ns=fwhm/bins_per_ns;
fprintf(1,'The center is at bin number %6.1f\n',xmax);
fprintf(1,'The FWHM is %6.1f picosec.\n',fwhm_ns*1000);