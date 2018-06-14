function centers = cntr3dg(img,pks,sz,idx2)

% fitting of 3D gaussian using custom algorithm based on 
% derivative-successive-approximation-least-squares.

% This fits the whole 3D gaussian, without summing rows or columns.

% Set up X and Y variables
% Be careful not to mix up rows/columns (need to check later)
imgsz=size(img);
X=1:imgsz(2); % let X be rows (need to check later)
Y=1:imgsz(1); % let Y be columns

x0=pks(1); % center of gaussian from "peaks"
y0=pks(2);
w=2;

% find array indices within Region of Interest (ROI)
XROIidx=((x0-sz)<=X) & (X<=(x0+sz));
YROIidx=((y0-sz)<=Y) & (Y<=(y0+sz));

% use only points within ROI for later calculations
XROI=X(XROIidx);
YROI=Y(YROIidx);
%XimgROI=img(circidx);
imgROI=img(YROI,XROI);

imgv=mat2col(imgROI); % flatten to column vector so can use mldivide

% make vectors:

x=XROI;
y=YROI;
[xx,yy]=meshgrid(x,y);
[xxv,xxdim]=mat2col(xx);
[yyv,yydim]=mat2col(yy); 

for cnt=1:50,
  
  xxvx0 = (xxv-x0);
  xxvx0_2 = xxvx0.^2;
  yyvy0 = (yyv-y0);
  yyvy0_2 = yyvy0.^2;
  x_2_y_2 = xxvx0_2 + yyvy0_2;
  
  g=exp(-x_2_y_2/2/w/w);
  gww=g/w/w;
  gpx=(xxv-x0).*gww;
  gpy=(yyv-y0).*gww;
  gpw=x_2_y_2.*gww/w;
  
  
%   g=exp(-((xxv-x0).^2+(yyv-y0).^2)/2/w/w);
%   gpx=(xxv-x0).*exp(-((xxv-x0).^2+(yyv-y0).^2)/2/w/w)/w/w;
%   gpy=(yyv-y0).*exp(-((xxv-x0).^2+(yyv-y0).^2)/2/w/w)/w/w;
%   gpw=((yyv-y0).^2+(xxv-x0).^2).*exp(-((xxv-x0).^2+(yyv-y0).^2)/2/w/w)/w/w/w;

        A=ones([length(g) 5]);

  A(:,1)=g;
  A(:,2)=gpx;
  A(:,3)=gpy;
  A(:,4)=gpw;
% www = warning ('on','all');
  c=A\imgv;
% warning(www);
% warning('on')
  dx=c(2)/c(1);
  dy=c(3)/c(1);
  dw=c(4)/c(1);
  ifit=c(1);

  % try not to overshoot
  if abs(dx)>.1;
    dx=dx/3;
  end
  if abs(dy)>.1;
    dy=dy/3;
  end
  if abs(dw)>.1;
    dw=dw/3;
  end
  
  % try not to WAY overshoot
  if abs(dx)>0.5,
    dx=0.2*sign(dx);
  end
  if abs(dy)>0.5,
    dy=0.2*sign(dy);
  end  

  x0=x0+dx;
  y0=y0+dy;
  w=w+dw;
  if (abs(dx)<1e-5) && (abs(dy)<1e-5) 
    break;
  end
  
    
end

% return a sane result: if position is >3 pixels from initial guess, 
% or is NaN, return initial guess
if abs(x0-pks(1))>3 || isnan(x0),
  x0=pks(1);
end
if abs(y0-pks(2))>3 || isnan(y0),
  y0=pks(2);
end

if w<1 || w>8 || isnan(w),
  w=2;
end

% recalculate the intensity, and integrate it:
g=exp(-((xxv-x0).^2+(yyv-y0).^2)/2/w/w);
A=ones([length(g) 2]);
A(:,1)=g;
c=A\imgv;
ifit=c(1)*sum(sum(g));
if isnan(ifit),
  ifit=sum(sum(imgv));
end


centers=[x0 y0 w^2 ifit];