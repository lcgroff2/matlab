
function H=mkffilt(img,bp1,bp2)

% MKFFILT.M

% function for making an ellipsoidal mask for a Fourier bandpass filter.

% Ellipsoidal instead of circular, because image is not usually square.

% img is the image file (used only to properly size the mask)

% bp1 is the low freq cutoff.  1 or 2 pixels, typical.  Non-integer
% values like 1.5 are OK.

% bp2 is the high frequency cutoff.  1 or 2 pixels typical.

% The output "H" is for passing to "FPASS" to do the actual
% filtering.

% H should look like an elliptical donut.

sz=size(img);

x=1:sz(1);
y=1:sz(2);

% center x,y
x=x-mean(x);
y=y-mean(y);

%xrad=-x(1)-bp2;
xrad=bp2;
%yrad=-y(1)-bp2;
yrad=bp2;

[xx,yy]=meshgrid(x,y);

% "normalized" radius for ellipse
Rsq=xx.*xx/xrad^2+yy.*yy/yrad^2;

H1=(1+exp(-30*(1-Rsq))).^(-1);

% regular radius in pixels

Rsq=xx.*xx+yy.*yy;

H2=exp(-Rsq/2/bp1);

H=H1-H2;

%imagesc(x,y,H)